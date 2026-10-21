classdef PhaseStability < handle
    % PHASESTABILITY General two-sided Michelsen tangent-plane-distance test
    % for an arbitrary feed (vapour-like, liquid-like or near-critical), FWI-only.
    %
    % Unlike solvers.StabilityTester (dew side only: vapour feed vs incipient
    % liquid), this tests the feed against BOTH vapour-like and liquid-like
    % trial phases, which a PT flash and a bubble-point seed require.
    %
    % Formulation (trial mole numbers W, single pressure P):
    %   h_i  = ln z_i + ln phi_i(z, P)
    %   SS   : ln W_i <- h_i - ln phi_i(x, P),       x = W / sum(W)
    %   tpd  = sum_i x_i [ln x_i + ln phi_i(x) - h_i]   (= -ln sum W at stationarity)
    % Unstable <=> tpd < -TPDTolerance at any non-trivial stationary point.
    %
    % Root selection: minimum-Gibbs root for feed and trial (required for a
    % valid tangent-plane criterion; same convention as StabilityTester).
    %
    % Confinement note: Wilson K-values carry no confinement information.
    % The default seed set is therefore widened (pure-heaviest and
    % pure-lightest trials) and callers should pass K-values from a nearby
    % converged confined state through 'ExtraK' whenever available; the
    % window between the bulk and confined saturation pressure is exactly
    % where a Wilson-only test is most likely to miss the second phase.

    properties (SetAccess = private)
        EOS thermo.ConfinedEOS
        MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 2000
        Tolerance     (1,1) double {mustBePositive} = 1e-10   % max |d ln W|
        TPDTolerance  (1,1) double {mustBePositive} = 1e-8
        TrivialTol    (1,1) double {mustBePositive} = 1e-4    % max |x - z| flagged as trivial
    end

    methods
        function obj = PhaseStability(eosEngine, opts)
            arguments
                eosEngine (1,1) thermo.ConfinedEOS
                opts.MaxIterations (1,1) double
                opts.Tolerance     (1,1) double
                opts.TPDTolerance  (1,1) double
                opts.TrivialTol    (1,1) double
            end
            obj.EOS = eosEngine;
            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
        end

        function res = test(obj, P, T, z, r_cap, opts)
            % P [Pa], T [K], z feed, r_cap [m] (Inf = bulk).
            % opts.ExtraK : [nc x m] K-values (y/x) used as additional seeds
            %               (both z.*K and z./K are tried for each column).
            % res.isUnstable, res.tpd, res.w (trial composition),
            % res.trialType ("liquid" | "vapor" | "none"),
            % res.K (y/x oriented initial K for a flash; NaN if stable),
            % res.seeds (per-seed diagnostics).
            arguments
                obj
                P     (1,1) double {mustBePositive}
                T     (1,1) double {mustBePositive}
                z     (:,1) double {mustBeNonnegative}
                r_cap (1,1) double {mustBePositive}
                opts.ExtraK double = []
            end
            z  = max(z, 1e-20);
            z  = z / sum(z);
            nc = numel(z);

            [lnphiZ, ~, ~] = obj.EOS.calculateState(P, T, z, 0, r_cap);
            h = log(z) + lnphiZ;

            seeds = obj.seedSet(P, T, z, opts.ExtraK);
            best  = struct('tpd', Inf, 'x', z, 'seed', 0);
            diag_ = repmat(struct('seed', 0, 'iterations', NaN, 'converged', false, ...
                                  'trivial', false, 'tpd', NaN), numel(seeds), 1);
            for s = 1:numel(seeds)
                try
                    r = obj.runSS(seeds{s}, z, h, P, T, r_cap);
                catch ME
                    if ~startsWith(string(ME.identifier), "ConfinedEOS:"), rethrow(ME); end
                    diag_(s).seed = s;
                    continue;
                end
                diag_(s) = struct('seed', s, 'iterations', r.iterations, 'converged', r.converged, ...
                                  'trivial', r.trivial, 'tpd', r.tpd);
                if ~r.trivial && r.tpd < best.tpd
                    best = struct('tpd', r.tpd, 'x', r.x, 'seed', s);
                end
            end

            res = struct('isUnstable', false, 'tpd', 0.0, 'w', z, 'trialType', "none", ...
                         'K', nan(nc, 1), 'bestSeed', 0, 'seeds', diag_);
            if ~isfinite(best.tpd)
                return;                                 % every seed collapsed onto the feed
            end
            res.tpd = best.tpd;
            res.w = best.x;
            res.bestSeed = best.seed;
            res.isUnstable = best.tpd < -obj.TPDTolerance;

            % Orientation: the trial heavier (by molar mass) than the feed is
            % the liquid-like incipient phase, so K = y/x = z/w; otherwise w is
            % vapour-like and K = w/z.
            MW = obj.EOS.Fluid.MW(:);
            if best.x.' * MW > z.' * MW
                res.trialType = "liquid";
                Kraw = z ./ max(best.x, 1e-300);
            else
                res.trialType = "vapor";
                Kraw = best.x ./ z;
            end
            if res.isUnstable
                res.K = Kraw;
            end
        end
    end

    methods (Access = private)
        function seeds = seedSet(obj, P, T, z, extraK)
            Kw = obj.EOS.wilsonK(P, T);
            [~, iH] = max(obj.EOS.Fluid.Tc);
            [~, iL] = min(obj.EOS.Fluid.Tc);
            heavy = 1e-3 * z;  heavy(iH) = 1;
            light = 1e-3 * z;  light(iL) = 1;
            seeds = {z .* Kw, z ./ Kw, heavy, light};
            if ~isempty(extraK)
                if size(extraK, 1) ~= numel(z)
                    extraK = extraK.';
                end
                for j = 1:size(extraK, 2)
                    Kj = extraK(:, j);
                    if all(isfinite(Kj)) && all(Kj > 0)
                        seeds(end+1:end+2) = {z .* Kj, z ./ Kj};
                    end
                end
            end
        end

        function r = runSS(obj, W0, z, h, P, T, r_cap)
            lnW = log(max(W0(:), 1e-300));
            converged = false;
            trivial   = false;
            it = 0;
            for it = 1:obj.MaxIterations
                x = exp(lnW - max(lnW));
                x = x / sum(x);
                [lnphi, ~, ~] = obj.EOS.calculateState(P, T, x, 0, r_cap);
                lnW_new = h - lnphi;
                err = max(abs(lnW_new - lnW));
                lnW = lnW_new;

                x = exp(lnW - max(lnW));
                x = x / sum(x);
                if max(abs(x - z)) < obj.TrivialTol
                    trivial = true;
                    break;
                end
                if err < obj.Tolerance
                    converged = true;
                    break;
                end
            end
            [lnphi, ~, ~] = obj.EOS.calculateState(P, T, x, 0, r_cap);
            tpd = x.' * (log(max(x, 1e-300)) + lnphi - h);
            if trivial
                tpd = 0.0;
            end
            r = struct('x', x, 'tpd', tpd, 'iterations', it, ...
                       'converged', converged, 'trivial', trivial);
        end
    end
end
