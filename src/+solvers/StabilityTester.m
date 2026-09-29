classdef StabilityTester < handle
    % STABILITYTESTER Tangent-plane-distance stability test of a vapor feed
    % against an incipient liquid, isobaric (FWI-only) or non-isobaric
    % (Combined: FWI + capillary pressure, P_L = P_V - Pcap).
    %
    % Formulation (Michelsen, trial mole numbers W):
    %   h_i  = ln z_i + ln phi_i^V(z, P_V) + ln P_V
    %   d_i  = ln W_i + ln phi_i^L(x, P_L) + ln P_L - h_i,   x = W / sum(W)
    %   tm   = 1 + sum_i W_i (d_i - 1)
    %   SS   : ln W_i <- h_i - ln phi_i^L(x, P_L) - ln P_L
    % At a stationary point d_i = 0 and tm = 1 - sum(W); instability <=> sum(W) > 1.
    % The reported tpd_min is the mole-fraction TPD, sum_i x_i d_i (= -ln sum W at
    % stationarity), identical in meaning to the previous implementation.
    %
    % In Combined mode Pcap(x) is updated with under-relaxation inside the SS
    % loop (lagged), using the shared ConfinedEOS.capillaryPressure.
    %
    % Root selection: minimum-Gibbs root for both feed and trial phase, as
    % required for a valid tangent-plane criterion. Forcing the vapor root on
    % the feed can report spurious instability; forcing a root on the trial
    % phase can miss a true one.

    properties (SetAccess = private)
        EOS thermo.ConfinedEOS
        MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 500
        Tolerance     (1,1) double {mustBePositive} = 1e-9    % max |d ln W| and |dPcap|/P_V
        TPDTolerance  (1,1) double {mustBePositive} = 1e-7    % unstable if tpd < -TPDTolerance
        TrivialTol    (1,1) double {mustBePositive} = 1e-4    % max |x - z| flagged as trivial
        RelaxPcap     (1,1) double {mustBeInRange(RelaxPcap, 0, 1, "exclude-lower")} = 0.5
        PressureFloor (1,1) double {mustBePositive} = 1e3     % [Pa] lower bound on P_L
    end

    methods
        function obj = StabilityTester(eosEngine, opts)
            arguments
                eosEngine (1,1) thermo.ConfinedEOS
                opts.MaxIterations (1,1) double
                opts.Tolerance     (1,1) double
                opts.TPDTolerance  (1,1) double
                opts.TrivialTol    (1,1) double
                opts.RelaxPcap     (1,1) double
                opts.PressureFloor (1,1) double
            end
            obj.EOS = eosEngine;
            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
        end

        function [isUnstable, w_star, tpd_min, K_factors, PV, PL, Pcap, info] = ...
                evaluateDewStability(obj, Pext, T, z, r_cap, opts)
            % Pext = P_V [Pa]; r_cap [m]; CapMode: "Combined" (legacy "Pc") or
            % "FWIOnly" (legacy "PL" / "NoCap").
            arguments
                obj
                Pext  (1,1) double {mustBePositive}
                T     (1,1) double {mustBePositive}
                z     (:,1) double {mustBeNonnegative}
                r_cap (1,1) double {mustBePositive}
                opts.CapMode (1,1) string = "Combined"
            end
            isCap = solvers.StabilityTester.isCapillaryActive(opts.CapMode, r_cap);

            z  = max(z, 1e-20);
            z  = z / sum(z);
            PV = Pext;
            nc = numel(z);

            % Feed (reference) state: minimum-Gibbs root
            [lnphiV, ~, V_V] = obj.EOS.calculateState(PV, T, z, 0, r_cap);
            h = log(z) + lnphiV + log(PV);

            % Trial seeds: Wilson liquid-like, and heaviest-component-rich
            Kw = obj.EOS.wilsonK(PV, T);
            [~, iHeavy] = max(obj.EOS.Fluid.Tc);
            seedHeavy = 1e-3 * z;
            seedHeavy(iHeavy) = 1;
            seeds = {z ./ Kw, seedHeavy};

            best = struct('tpd', Inf);
            info = struct('seed', {}, 'iterations', {}, 'converged', {}, ...
                          'trivial', {}, 'tpd', {}, 'tm', {});
            for s = 1:numel(seeds)
                try
                    r = obj.runSS(seeds{s}, z, h, T, PV, V_V, r_cap, isCap, nc);
                catch ME
                    if ~startsWith(ME.identifier, "ConfinedEOS:"), rethrow(ME); end
                    info(s) = struct('seed', s, 'iterations', NaN, 'converged', false, ...
                                     'trivial', false, 'tpd', NaN, 'tm', NaN);
                    continue;
                end
                info(s) = struct('seed', s, 'iterations', r.iterations, 'converged', r.converged, ...
                                 'trivial', r.trivial, 'tpd', r.tpd, 'tm', r.tm);
                if ~r.trivial && r.tpd < best.tpd
                    best = r;
                end
            end

            if ~isfinite(best.tpd)
                % Every seed collapsed onto the feed: stable
                isUnstable = false;
                w_star     = z;
                tpd_min    = 0.0;
                PL         = PV;
                Pcap       = 0.0;
                K_factors  = nan(nc, 1);
                return;
            end

            tpd_min    = best.tpd;
            w_star     = best.x;
            PL         = best.PL;
            Pcap       = best.Pcap;
            isUnstable = tpd_min < -obj.TPDTolerance;
            if isUnstable
                K_factors = z ./ max(w_star, 1e-14);      % K = y/x (dew tracing)
            else
                K_factors = nan(nc, 1);
            end
        end
    end

    methods (Static)
        function tf = isCapillaryActive(capMode, r_cap)
            % Single interpretation of the capillary-mode switch for all solvers.
            m = lower(string(capMode));
            if any(m == ["combined", "pc"])
                tf = isfinite(r_cap);
            elseif any(m == ["fwionly", "pl", "nocap"])
                tf = false;
            else
                error('StabilityTester:UnknownCapMode', ...
                    'CapMode "%s" not recognised (use "Combined" or "FWIOnly").', capMode);
            end
        end
    end

    methods (Access = private)
        function r = runSS(obj, W0, z, h, T, PV, V_V, r_cap, isCap, nc)
            lnW  = log(max(W0(:), 1e-300));
            Pcap = 0.0;
            PL   = PV;
            converged = false;
            trivial   = false;

            for it = 1:obj.MaxIterations
                x = exp(lnW - max(lnW));
                x = x / sum(x);

                [lnphiL, ~, V_L] = obj.EOS.calculateState(PL, T, x, 0, r_cap);
                lnW_new = h - lnphiL - log(PL);

                % Lagged capillary update (relaxed)
                dPcap = 0.0;
                if isCap
                    if max(abs(x - z)) < obj.TrivialTol
                        Pcap_new = 0.0;             % no interface for a trivial trial phase
                    else
                        Pcap_new = obj.EOS.capillaryPressure(x, z, V_L, V_V, r_cap);
                    end
                    Pcap_upd = (1 - obj.RelaxPcap) * Pcap + obj.RelaxPcap * Pcap_new;
                    dPcap    = abs(Pcap_upd - Pcap) / PV;
                    Pcap     = Pcap_upd;
                    PL       = max(PV - Pcap, obj.PressureFloor);
                end

                err = max([max(abs(lnW_new - lnW)), dPcap]);
                lnW = lnW_new;

                x_new = exp(lnW - max(lnW));
                x_new = x_new / sum(x_new);
                if max(abs(x_new - z)) < obj.TrivialTol
                    trivial = true;
                    break;
                end
                if err < obj.Tolerance
                    converged = true;
                    break;
                end
            end

            % Evaluate tm and tpd consistently at the final W and P_L
            W = exp(lnW);
            x = W / sum(W);
            [lnphiL, ~, ~] = obj.EOS.calculateState(PL, T, x, 0, r_cap);
            d  = log(W) + lnphiL + log(PL) - h;
            tm = 1 + sum(W .* (d - 1));
            dx = log(x) + lnphiL + log(PL) - h;         % mole-fraction TPD
            tpd = x.' * dx;

            if trivial
                tpd = 0.0;
            end
            r = struct('x', x, 'tpd', tpd, 'tm', tm, 'PL', PL, 'Pcap', Pcap, ...
                       'iterations', it, 'converged', converged, 'trivial', trivial, ...
                       'nc', nc);
        end
    end
end
