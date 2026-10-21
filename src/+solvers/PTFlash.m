classdef PTFlash < handle
    % PTFLASH Isothermal two-phase (VL) PT flash for confined mixtures, FWI-only.
    %
    % Sequence: PhaseStability (two-sided TPD) -> successive substitution on
    % ln K with a bounded Rachford-Rice solve -> Newton on ln K (central-
    % difference Jacobian, Armijo line search) -> acceptance checks.
    %
    %   Unknowns   ln K_i,      K_i = y_i / x_i
    %   Residual   F_i = ln K_i - [ln phi_i^L(x, P) - ln phi_i^V(y, P)]
    %   with x, y from Rachford-Rice: sum_i z_i (K_i - 1) / (1 + beta (K_i - 1)) = 0
    %
    % beta is the vapour MOLE fraction. Root selection is minimum-Gibbs for
    % both phases (Michelsen); phases are labelled by mass density.
    %
    % Scope: FWI-only (one pressure for both phases). The Combined (Pcap)
    % mode needs a reference-pressure convention for a single-phase liquid
    % feed (gas or liquid pressure) before a PT flash is well defined; it is
    % intentionally not implemented here (see the dew-side FlashEngine).
    %
    % Results are accepted only if the two-phase Gibbs energy is below the
    % single-phase value, 0 < beta < 1, and the split is non-trivial.

    properties (SetAccess = private)
        EOS       thermo.ConfinedEOS
        Stability solvers.PhaseStability

        SSMaxIterations (1,1) double {mustBeInteger, mustBePositive} = 300
        SwitchTol       (1,1) double {mustBePositive} = 1e-5    % SS -> Newton when max|d ln K| < SwitchTol
        MaxIterations   (1,1) double {mustBeInteger, mustBePositive} = 60
        Tolerance       (1,1) double {mustBePositive} = 1e-11   % ||F||_inf
        JacEpsilon      (1,1) double {mustBePositive} = 1e-6
        MaxStepLnK      (1,1) double {mustBePositive} = 0.5
        LineSearchC1    (1,1) double {mustBePositive} = 1e-4
        MaxLineSearch   (1,1) double {mustBeInteger, mustBePositive} = 20
        TrivialTol      (1,1) double {mustBePositive} = 1e-4    % ||ln K||_2 below which the split is trivial
        GibbsTol        (1,1) double {mustBePositive} = 1e-12
    end

    methods
        function obj = PTFlash(eosEngine, opts)
            arguments
                eosEngine (1,1) thermo.ConfinedEOS
                opts.Stability solvers.PhaseStability = solvers.PhaseStability.empty
                opts.SSMaxIterations (1,1) double
                opts.SwitchTol       (1,1) double
                opts.MaxIterations   (1,1) double
                opts.Tolerance       (1,1) double
                opts.JacEpsilon      (1,1) double
                opts.MaxStepLnK      (1,1) double
                opts.LineSearchC1    (1,1) double
                opts.MaxLineSearch   (1,1) double
                opts.TrivialTol      (1,1) double
                opts.GibbsTol        (1,1) double
            end
            obj.EOS = eosEngine;
            if isempty(opts.Stability)
                obj.Stability = solvers.PhaseStability(eosEngine);
            else
                if opts.Stability.EOS ~= eosEngine
                    error('PTFlash:EngineMismatch', 'PhaseStability must wrap the same ConfinedEOS instance.');
                end
                obj.Stability = opts.Stability;
            end
            opts = rmfield(opts, 'Stability');
            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
        end

        function res = solve(obj, P, T, z, r_cap, opts)
            % P [Pa], T [K], z feed, r_cap [m] (Inf = bulk).
            % opts.K_seed : K = y/x from a nearby converged state (warm start;
            %               also passed to the stability test as an extra seed).
            % opts.SkipStability : trust K_seed and go straight to the split
            %               (only for marching along a known two-phase path).
            arguments
                obj
                P     (1,1) double {mustBePositive}
                T     (1,1) double {mustBePositive}
                z     (:,1) double {mustBeNonnegative}
                r_cap (1,1) double {mustBePositive}
                opts.K_seed double = []
                opts.SkipStability (1,1) logical = false
            end
            nc = obj.EOS.Fluid.NC;
            if numel(z) ~= nc
                error('PTFlash:CompositionSize', 'z has %d entries; NC = %d.', numel(z), nc);
            end
            z = max(z, 1e-20);
            z = z / sum(z);
            Kseed = opts.K_seed(:);

            % ---------- 1. Stability ----------
            st = struct('isUnstable', true, 'tpd', NaN, 'trialType', "seed");
            if opts.SkipStability && ~isempty(Kseed)
                K0 = Kseed;
            else
                st = obj.Stability.test(P, T, z, r_cap, 'ExtraK', Kseed);
                if ~st.isUnstable
                    res = obj.singlePhase(P, T, z, r_cap, "stable", st);
                    return;
                end
                K0 = st.K;
            end

            % ---------- 2. Successive substitution ----------
            lnK = log(max(K0, 1e-300));
            reason = "ss-max-iterations";
            itSS = 0;
            try
                for itSS = 1:obj.SSMaxIterations
                    [lnKnew, ~] = obj.update(lnK, P, T, z, r_cap);
                    err = max(abs(lnKnew - lnK));
                    lnK = lnKnew;
                    if norm(lnK) < obj.TrivialTol
                        reason = "trivial";
                        break;
                    end
                    if err < obj.SwitchTol
                        reason = "switch";
                        break;
                    end
                end
            catch ME
                if ~obj.isRecoverable(ME), rethrow(ME); end
                reason = "model-failure-ss";
            end
            if reason == "trivial" || reason == "model-failure-ss"
                res = obj.singlePhase(P, T, z, r_cap, reason, st);
                return;
            end

            % ---------- 3. Newton on ln K ----------
            [lnK, nIt, nReason] = obj.newton(lnK, P, T, z, r_cap);

            % ---------- 4. Acceptance ----------
            res = obj.package(lnK, P, T, z, r_cap, st);
            res.iterationsSS = itSS;
            res.iterationsNewton = nIt;
            res.reason = nReason;
            res.converged = nReason == "converged";
            if ~res.converged
                return;
            end
            if norm(lnK) < obj.TrivialTol
                res = obj.singlePhase(P, T, z, r_cap, "trivial", st);
            elseif res.beta <= 0 || res.beta >= 1
                res = obj.singlePhase(P, T, z, r_cap, "negative-flash", st);
            elseif res.dG >= -obj.GibbsTol
                res = obj.singlePhase(P, T, z, r_cap, "gibbs-rejected", st);
            end
        end
    end

    methods (Access = private)
        function [lnKnew, s] = update(obj, lnK, P, T, z, r_cap)
            K = exp(lnK);
            beta = solvers.PTFlash.rachfordRice(z, K);
            [x, y] = solvers.PTFlash.phaseCompositions(z, K, beta);
            [lnphiL, ZL, VL] = obj.EOS.calculateState(P, T, x, 0, r_cap);
            [lnphiV, ZV, VV] = obj.EOS.calculateState(P, T, y, 0, r_cap);
            lnKnew = lnphiL - lnphiV;
            s = struct('beta', beta, 'x', x, 'y', y, 'lnphiL', lnphiL, 'lnphiV', lnphiV, ...
                       'ZL', ZL, 'ZV', ZV, 'VL', VL, 'VV', VV);
        end

        function [F, ok] = residual(obj, lnK, P, T, z, r_cap)
            try
                [lnKnew, ~] = obj.update(lnK, P, T, z, r_cap);
                F = lnK - lnKnew;
                ok = all(isfinite(F));
            catch ME
                if ~obj.isRecoverable(ME), rethrow(ME); end
                F = [];  ok = false;
            end
        end

        function [lnK, it, reason] = newton(obj, lnK, P, T, z, r_cap)
            [F, ok] = obj.residual(lnK, P, T, z, r_cap);
            reason = "max-iterations";
            it = 0;
            if ~ok
                reason = "initial-state-failure";
                return;
            end
            n = numel(lnK);
            meritHist = 0.5 * (F.' * F);
            for it = 1:obj.MaxIterations
                if norm(F, Inf) < obj.Tolerance
                    reason = "converged";
                    return;
                end
                J = zeros(n);
                for j = 1:n
                    h  = obj.JacEpsilon * max(1, abs(lnK(j)));
                    up = lnK; up(j) = up(j) + h;
                    um = lnK; um(j) = um(j) - h;
                    [Fp, okp] = obj.residual(up, P, T, z, r_cap);
                    [Fm, okm] = obj.residual(um, P, T, z, r_cap);
                    if okp && okm
                        J(:, j) = (Fp - Fm) / (2 * h);
                    elseif okp
                        J(:, j) = (Fp - F) / h;
                    elseif okm
                        J(:, j) = (F - Fm) / h;
                    end
                end
                if rcond(J) > 1e-14
                    d = -(J \ F);
                else
                    d = -((J.' * J + 1e-10 * max(1, norm(J, 1)^2) * eye(n)) \ (J.' * F));
                end
                sMax = norm(d, Inf);
                if sMax > obj.MaxStepLnK
                    d = d * (obj.MaxStepLnK / sMax);
                end
                % Non-monotone reference with absolute floor (near-critical splits)
                fRef = max(max(meritHist(max(1, end - 4):end)), 1e-8);
                alpha = 1.0;
                accepted = false;
                for ls = 1:obj.MaxLineSearch
                    [Ft, okt] = obj.residual(lnK + alpha * d, P, T, z, r_cap);
                    if okt && 0.5 * (Ft.' * Ft) <= (1 - 2 * obj.LineSearchC1 * alpha) * fRef
                        accepted = true;
                        break;
                    end
                    alpha = 0.5 * alpha;
                end
                if ~accepted
                    reason = "line-search-failure";
                    return;
                end
                lnK = lnK + alpha * d;
                F = Ft;
                meritHist(end+1) = 0.5 * (F.' * F); %#ok<AGROW>
            end
            if norm(F, Inf) < obj.Tolerance
                reason = "converged";
            end
        end

        function res = package(obj, lnK, P, T, z, r_cap, st)
            [~, s] = obj.update(lnK, P, T, z, r_cap);
            MW  = obj.EOS.Fluid.MW(:);
            [lnphiZ, ~, ~] = obj.EOS.calculateState(P, T, z, 0, r_cap);
            gL = s.x.' * (log(max(s.x, 1e-300)) + s.lnphiL);
            gV = s.y.' * (log(max(s.y, 1e-300)) + s.lnphiV);
            gZ = z.' * (log(z) + lnphiZ);
            dG = (1 - s.beta) * gL + s.beta * gV - gZ;          % dimensionless G/(RT) change

            rhoL = 1e-3 * (s.x.' * MW) / s.VL;                  % [kg/m^3]
            rhoV = 1e-3 * (s.y.' * MW) / s.VV;
            % Phase identity by MASS density (consistent with FlashEngine)
            if rhoV > rhoL
                [s.x, s.y] = deal(s.y, s.x);
                [s.ZL, s.ZV] = deal(s.ZV, s.ZL);
                [s.VL, s.VV] = deal(s.VV, s.VL);
                [rhoL, rhoV] = deal(rhoV, rhoL);
                s.beta = 1 - s.beta;
                lnK = -lnK;
            end
            volL = (1 - s.beta) * s.VL;
            volV = s.beta * s.VV;
            res = struct('phases', 2, 'converged', true, 'reason', "", ...
                'P', P, 'T', T, 'r_cap', r_cap, 'z', z, ...
                'beta', s.beta, 'x', s.x, 'y', s.y, 'K', exp(lnK), ...
                'Z_L', s.ZL, 'Z_V', s.ZV, 'V_L', s.VL, 'V_V', s.VV, ...
                'rhoM_L', rhoL, 'rhoM_V', rhoV, ...
                'V_total', volL + volV, ...                      % [m^3 per mol feed]
                'SL', volL / (volL + volV), ...                  % liquid volume fraction
                'dG', dG, 'tpd', st.tpd, 'trialType', st.trialType, ...
                'label', "two-phase", 'iterationsSS', 0, 'iterationsNewton', 0);
        end

        function res = singlePhase(obj, P, T, z, r_cap, reason, st)
            [~, Z, V] = obj.EOS.calculateState(P, T, z, 0, r_cap);
            MW = obj.EOS.Fluid.MW(:);
            label = obj.singlePhaseLabel(T, z);
            isV = label == "vapor";
            nc = numel(z);
            res = struct('phases', 1, 'converged', true, 'reason', reason, ...
                'P', P, 'T', T, 'r_cap', r_cap, 'z', z, ...
                'beta', double(isV), 'x', z, 'y', z, 'K', ones(nc, 1), ...
                'Z_L', Z, 'Z_V', Z, 'V_L', V, 'V_V', V, ...
                'rhoM_L', 1e-3 * (z.' * MW) / V, 'rhoM_V', 1e-3 * (z.' * MW) / V, ...
                'V_total', V, 'SL', double(~isV), ...
                'dG', 0.0, 'tpd', st.tpd, 'trialType', st.trialType, ...
                'label', label, 'iterationsSS', 0, 'iterationsNewton', 0);
        end

        function label = singlePhaseLabel(obj, T, z)
            % Li's pseudo-critical temperature with BULK critical properties.
            % A label only (relative permeability bookkeeping); it carries no
            % confinement information, so do not use it near the confined
            % critical region without checking the phase envelope.
            Vc = obj.EOS.Fluid.Vc(:);
            Tc = obj.EOS.Fluid.Tc(:);
            Tpc = (z .* Vc).' * Tc / (z.' * Vc);
            if T > Tpc
                label = "vapor";
            else
                label = "liquid";
            end
        end
    end

    methods (Static)
        function beta = rachfordRice(z, K)
            % Vapour mole fraction from Rachford-Rice, safeguarded Newton-
            % bisection on the monotone function g(beta) inside the asymptote
            % window (1/(1-Kmax), 1/(1-Kmin)). Negative-flash values outside
            % [0, 1] are returned unclipped; the caller decides.
            Kmax = max(K);  Kmin = min(K);
            if Kmax <= 1
                beta = -Inf;  return;            % all K <= 1: liquid side
            elseif Kmin >= 1
                beta = Inf;   return;            % all K >= 1: vapour side
            end
            lo = 1 / (1 - Kmax);
            hi = 1 / (1 - Kmin);
            Km1 = K - 1;
            beta = 0.5;
            if ~(beta > lo && beta < hi), beta = 0.5 * (lo + hi); end
            for it = 1:200
                den = 1 + beta * Km1;
                g  = sum(z .* Km1 ./ den);
                dg = -sum(z .* Km1.^2 ./ den.^2);
                if g > 0, lo = beta; else, hi = beta; end
                bn = beta - g / dg;
                if ~(bn > lo && bn < hi)
                    bn = 0.5 * (lo + hi);
                end
                if abs(bn - beta) <= 1e-15 * max(1, abs(beta))
                    beta = bn;
                    return;
                end
                beta = bn;
            end
        end

        function [x, y] = phaseCompositions(z, K, beta)
            if ~isfinite(beta)
                x = z;  y = z;  return;
            end
            x = z ./ (1 + beta * (K - 1));
            y = K .* x;
            x = x / sum(x);
            y = y / sum(y);
        end
    end

    methods (Static, Access = private)
        function tf = isRecoverable(ME)
            tf = startsWith(string(ME.identifier), "ConfinedEOS:");
        end
    end
end
