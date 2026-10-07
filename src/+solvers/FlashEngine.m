classdef FlashEngine < handle
    % FLASHENGINE Dew-point solver for confined mixtures (incipient-liquid
    % formulation), isobaric (FWI-only) or non-isobaric (Combined: FWI + Pcap).
    %
    % Unknowns   u = [ln K (nc); ln P_V; ln Pcap (Combined only)]
    % Residuals  F_i    = ln K_i - [ln phi_i^L(x, P_L) + ln P_L - ln phi_i^V(z, P_V) - ln P_V]
    %            F_nc+1 = sum_i z_i / K_i - 1
    %            F_nc+2 = Pcap_YL(x, z, V_L, V_V) / Pcap - 1           (Combined)
    % with x = (z./K) / sum(z./K), P_L = P_V - Pcap.
    %
    % Solvers: successive substitution (standalone or pre-conditioner),
    % Newton-Raphson and Broyden (good update on J), both globalised with an
    % Armijo backtracking line search on 0.5*||F||^2.
    %
    % Root selection: forced (+1 liquid, -1 vapor) roots, as required by the
    % incipient-phase formulation to keep the phases distinct. Trivial and
    % coalesced-root states are rejected, never deflected.

    properties (Constant)
        psi2Pa = 6894.757
    end

    properties (SetAccess = private)
        EOS       thermo.ConfinedEOS
        Stability solvers.StabilityTester

        % Newton / Broyden
        MaxIterations   (1,1) double {mustBeInteger, mustBePositive} = 200
        Tolerance       (1,1) double {mustBePositive} = 1e-8     % ||F||_inf
        JacEpsilon      (1,1) double {mustBePositive} = 1e-6
        LineSearchC1    (1,1) double {mustBePositive} = 1e-4
        MaxLineSearch   (1,1) double {mustBeInteger, mustBePositive} = 14
        JacRegularizer  (1,1) double {mustBePositive} = 1e-10    % LM damping when J is ill-conditioned
        MaxStepSizeK    (1,1) double {mustBePositive} = 0.20     % max |d ln K|
        MaxStepSizeP    (1,1) double {mustBePositive} = 0.05     % max |d ln P_V|
        MaxStepSizePcap (1,1) double {mustBePositive} = 0.50     % max |d ln Pcap|

        % Successive substitution
        PreconditionWithSS (1,1) logical = true
        SSMaxIterations    (1,1) double {mustBeInteger, mustBePositive} = 35
        SSTolerance        (1,1) double {mustBePositive} = 1e-5
        SSDamping          (1,1) double {mustBeInRange(SSDamping, 0, 1, "exclude-lower")} = 0.70
        SSMaxStepP         (1,1) double {mustBePositive} = log(1.1)
        DewBranch          (1,1) string {mustBeMember(DewBranch, ["upper","lower"])} = "upper"

        % Physical guards
        TrivialTol      (1,1) double {mustBePositive} = 0.01     % ||ln K||_2 below which a state is trivial
        SeedMinLnK      (1,1) double {mustBePositive} = 0.05     % minimum ||ln K||_2 to accept a TPD seed
        CoalescenceTol  (1,1) double {mustBePositive} = 1e-5     % |Z_V - Z_L|
        BackoffFactor   (1,1) double {mustBeInRange(BackoffFactor, 0, 1, "exclusive")} = 0.98
        MaxBackoffs     (1,1) double {mustBeInteger, mustBePositive} = 50
        PressureMin     (1,1) double {mustBePositive} = 1e5                  % [Pa]
        PressureMax     (1,1) double {mustBePositive} = 35000 * 6894.757     % [Pa]
        PressureFloor   (1,1) double {mustBePositive} = 1e3                  % [Pa] floor on P_L
        % Phase-identity test at the returned root. Never molar density or Z:
        % near the critical point of an asymmetric mixture the incipient liquid
        % can have the higher Z while being heavier and denser by MASS.
        %   "massDensity": rho_m,L > rho_m,V   (rho_m = M_mix / V)
        %   "composition": incipient phase enriched in the heaviest component
        %   "both" | "none"
        PhaseIdentityTest (1,1) string {mustBeMember(PhaseIdentityTest, ...
            ["massDensity", "composition", "both", "none"])} = "massDensity" 
    end

    methods
        function obj = FlashEngine(eosEngine, stabilityEngine, opts)
            arguments
                eosEngine       (1,1) thermo.ConfinedEOS
                stabilityEngine (1,1) solvers.StabilityTester
                opts.MaxIterations
                opts.Tolerance
                opts.JacEpsilon
                opts.LineSearchC1
                opts.MaxLineSearch
                opts.JacRegularizer
                opts.MaxStepSizeK
                opts.MaxStepSizeP
                opts.MaxStepSizePcap
                opts.PreconditionWithSS
                opts.SSMaxIterations
                opts.SSTolerance
                opts.SSDamping
                opts.SSMaxStepP
                opts.DewBranch
                opts.TrivialTol
                opts.SeedMinLnK
                opts.CoalescenceTol
                opts.BackoffFactor
                opts.MaxBackoffs
                opts.PhaseIdentityTest
                opts.PressureMin
                opts.PressureMax
                opts.PressureFloor
            end
            if stabilityEngine.EOS ~= eosEngine
                error('FlashEngine:EngineMismatch', ...
                    'StabilityTester must wrap the same ConfinedEOS instance.');
            end
            obj.EOS = eosEngine;
            obj.Stability = stabilityEngine;
            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
        end

        function [Pdew, K_final, Pcap_final, Pliq_final, solverStats] = ...
                solveDewPoint(obj, T, P_guess, z, r_cap, opts)
            arguments
                obj
                T       (1,1) double {mustBePositive}
                P_guess (1,1) double {mustBePositive}
                z       (:,1) double {mustBeNonnegative}
                r_cap   (1,1) double {mustBePositive}
                opts.Solver         (1,1) string = "newton"
                opts.CapMode        (1,1) string = "Combined"
                opts.K_seed         double = []
                opts.Pcap_seed      double = []
                opts.PreconditionSS (1,1) logical = obj.PreconditionWithSS
                opts.SSMaxIter      (1,1) double = obj.SSMaxIterations
                opts.SSTol          (1,1) double = obj.SSTolerance
            end
            solver = lower(opts.Solver);
            if solver == "broyden", solver = "quasinewton"; end
            if ~any(solver == ["newton", "quasinewton", "ss"])
                error('FlashEngine:UnknownSolver', 'Solver "%s" not recognised.', opts.Solver);
            end

            z  = z / sum(z);
            nc = obj.EOS.Fluid.NC;
            if numel(z) ~= nc
                error('FlashEngine:CompositionSize', 'z has %d entries; NC = %d.', numel(z), nc);
            end
            isCap = solvers.StabilityTester.isCapillaryActive(opts.CapMode, r_cap);
            if isCap && cosd(obj.EOS.Rock.Theta) <= 0
                error('FlashEngine:NonWettingLiquid', ...
                    'theta >= 90 deg gives Pcap <= 0; the ln(Pcap) formulation does not apply.');
            end

            % ---------- 1. Initial estimates ----------
            [K0, Pcap0, seedSource] = obj.initialEstimates(T, P_guess, z, r_cap, isCap, opts);
            [~, iHeavy] = max(obj.EOS.Fluid.Tc);
            [~, iLight] = min(obj.EOS.Fluid.Tc);
            K0(iHeavy) = min(K0(iHeavy), 0.95);
            K0(iLight) = max(K0(iLight), 1.02);
            PV0 = obj.clampP(P_guess);

            % ---------- 2. Successive substitution ----------
            ssStats = [];
            if solver == "ss" || opts.PreconditionSS
                if solver == "ss"
                    [maxIt, tol] = deal(obj.MaxIterations, obj.Tolerance);
                else
                    [maxIt, tol] = deal(opts.SSMaxIter, opts.SSTol);
                end
                [PV_ss, K_ss, Pcap_ss, ssStats] = ...
                    obj.executeSS(T, PV0, z, r_cap, K0, Pcap0, isCap, maxIt, tol);

                if solver == "ss"
                    Pdew        = PV_ss;
                    K_final     = K_ss;
                    Pcap_final  = Pcap_ss;
                    Pliq_final  = obj.liquidPressure(PV_ss, Pcap_ss, isCap);
                    solverStats = ssStats;
                    solverStats.seed = seedSource;
                    return;
                end
                if ssStats.usable
                    [PV0, K0, Pcap0] = deal(PV_ss, K_ss, Pcap_ss);
                end
            end

            % ---------- 3. Newton / Broyden ----------
            if isCap
                if Pcap0 <= 0
                    Pcap0 = obj.consistentPcap(T, PV0, z, r_cap, K0);
                end
                u0 = [log(max(K0, 1e-14)); log(PV0); log(max(Pcap0, 1.0))];
            else
                u0 = [log(max(K0, 1e-14)); log(PV0)];
            end
            [u, solverStats] = obj.executeNewton(u0, T, z, r_cap, isCap, solver == "quasinewton");

            K_final = exp(u(1:nc));
            Pdew    = exp(u(nc+1));
            if isCap
                Pcap_final = exp(u(nc+2));
            else
                Pcap_final = 0.0;
            end
            Pliq_final = obj.liquidPressure(Pdew, Pcap_final, isCap);
            solverStats.seed = seedSource;
            solverStats.ss   = ssStats;
        end
    end

    methods (Access = private)
        % ================= Initialisation =================
        function [K0, Pcap0, source] = initialEstimates(obj, T, P_guess, z, r_cap, isCap, opts)
            nc = numel(z);
            if ~isempty(opts.K_seed)
                if numel(opts.K_seed) ~= nc
                    error('FlashEngine:SeedSize', 'K_seed must have NC = %d entries.', nc);
                end
                K0 = opts.K_seed(:);
                Pcap0 = 0.0;
                if isCap && ~isempty(opts.Pcap_seed)
                    Pcap0 = max(opts.Pcap_seed, 0.0);
                end
                source = "user";
                return;
            end

            capMode = "FWIOnly";
            if isCap, capMode = "Combined"; end
            try
                [isUnstable, ~, ~, K_tpd, ~, ~, Pcap_tpd] = ...
                    obj.Stability.evaluateDewStability(P_guess, T, z, r_cap, 'CapMode', capMode);
                if isUnstable && all(isfinite(K_tpd)) && norm(log(K_tpd)) > obj.SeedMinLnK
                    K0 = K_tpd(:);
                    Pcap0 = isCap * Pcap_tpd;
                    source = "tpd";
                    return;
                end
            catch ME
                if ~obj.isModelFailure(ME), rethrow(ME); end
            end
            K0 = obj.EOS.wilsonK(P_guess, T);
            Pcap0 = 0.0;
            source = "wilson";
        end

        function Pcap = consistentPcap(obj, T, PV, z, r_cap, K)
            % Few fixed-point passes so ln(Pcap) starts on the right scale.
            x = z ./ K;  x = x / sum(x);
            Pcap = 0.0;
            try
                [~, ~, V_V] = obj.EOS.calculateState(PV, T, z, -1, r_cap);
                for k = 1:5
                    PL = obj.liquidPressure(PV, Pcap, true);
                    [~, ~, V_L] = obj.EOS.calculateState(PL, T, x, +1, r_cap);
                    Pcap = obj.EOS.capillaryPressure(x, z, V_L, V_V, r_cap);
                end
            catch ME
                if ~obj.isModelFailure(ME), rethrow(ME); end
            end
            Pcap = min(max(Pcap, 1.0), 0.5 * PV);
        end

        % ================= Successive substitution =================
        function [PV, K, Pcap, stats] = executeSS(obj, T, PV0, z, r_cap, K0, Pcap0, isCap, maxIter, tol)
            lnK  = log(max(K0(:), 1e-14));
            lnP  = log(obj.clampP(PV0));
            Pcap = isCap * max(Pcap0, 0.0);
            branchSlope = -1;                         % d ln S / d ln P_V on the upper dew branch
            if obj.DewBranch == "lower", branchSlope = +1; end

            prev = [];
            nBack = 0;
            converged = false;
            reason = "max-iterations";
            err = Inf;
            st = obj.phaseState(z, z, NaN, NaN, NaN, NaN, NaN);

            for iter = 1:maxIter
                PV = exp(lnP);
                PL = obj.liquidPressure(PV, Pcap, isCap);
                x  = z .* exp(-lnK);
                x  = x / sum(x);

                try
                    [lnphiV, ZV, VV] = obj.EOS.calculateState(PV, T, z, -1, r_cap);
                    [lnphiL, ZL, VL] = obj.EOS.calculateState(PL, T, x, +1, r_cap);
                catch ME
                    if ~obj.isModelFailure(ME), rethrow(ME); end
                    [lnP, nBack, prev] = obj.backoff(lnP, nBack);
                    if nBack > obj.MaxBackoffs, reason = "model-failure"; break; end
                    continue;
                end

                lnK_new = lnphiL + log(PL) - lnphiV - log(PV);
                if norm(lnK_new) < obj.TrivialTol || abs(ZV - ZL) < obj.CoalescenceTol
                    [lnP, nBack, prev] = obj.backoff(lnP, nBack);
                    if nBack > obj.MaxBackoffs, reason = "trivial"; break; end
                    continue;
                end
                st = obj.phaseState(z, x, ZV, ZL, VV, VL, sum(z .* exp(-lnK_new)) - 1);

                % Pressure update: secant on ln S(ln P_V) = 0, branch-signed start
                lnS = log(sum(z .* exp(-lnK_new)));
                slope = branchSlope;
                if ~isempty(prev) && abs(lnP - prev(1)) > 1e-12
                    sec = (lnS - prev(2)) / (lnP - prev(1));
                    if isfinite(sec) && abs(sec) > 1e-3
                        slope = sec;
                    end
                end
                dlnP = max(min(-lnS / slope, obj.SSMaxStepP), -obj.SSMaxStepP);
                prev = [lnP, lnS];

                % Capillary pressure from the current incipient-liquid state
                errCap = 0.0;
                if isCap
                    Pcap_new = obj.EOS.capillaryPressure(x, z, VL, VV, r_cap);
                    errCap   = abs(Pcap_new - Pcap) / max(Pcap, obj.PressureFloor);
                    Pcap     = obj.SSDamping * Pcap_new + (1 - obj.SSDamping) * Pcap;
                end

                err = max([norm(lnK_new - lnK, Inf), abs(lnS), errCap]);
                lnK = obj.SSDamping * lnK_new + (1 - obj.SSDamping) * lnK;
                lnP = log(obj.clampP(exp(lnP + dlnP)));

                if err < tol
                    converged = true;
                    reason = "converged";
                    break;
                end
            end

            PV = exp(lnP);
            K  = exp(lnK);
            usable = all(isfinite(K)) && isfinite(PV) && norm(lnK) >= obj.TrivialTol ...
                     && ~any(reason == ["model-failure", "trivial"]);
            stats = struct('converged', converged, 'usable', usable, 'reason', reason, ...
                'iterations', iter, 'backoffs', nBack, 'residual_norm', err, ...
                'solver', "SS-" + obj.modeLabel(isCap), ...
                'Z_V', st.Z_V, 'Z_L', st.Z_L, 'rho_v', 1 / st.V_V, 'rho_l', 1 / st.V_L, ...
                'sumX_res', st.sumX_res, 'rhoM_V', st.rhoM_V, 'rhoM_L', st.rhoM_L, ...
                'heavyEnriched', st.heavyEnriched, 'incipientIsLiquid', obj.identityOK(st), ...
                'Pcap', Pcap, 'PL', obj.liquidPressure(PV, Pcap, isCap));
        end

        function [lnP, nBack, prev] = backoff(obj, lnP, nBack)
            lnP   = log(obj.clampP(exp(lnP) * obj.BackoffFactor));
            nBack = nBack + 1;
            prev  = [];                               % secant history invalid after a jump
        end

        % ================= Newton / Broyden =================
        function [u, stats] = executeNewton(obj, u0, T, z, r_cap, isCap, useBroyden)
            nc = obj.EOS.Fluid.NC;
            u  = obj.project(u0, isCap);
            [F, st, ok] = obj.residual(u, T, z, r_cap, isCap);
            if ~ok
                stats = obj.newtonStats(false, "initial-state-failure", 0, Inf, st, useBroyden, u, nc, isCap);
                return;
            end

            J = [];
            freshJ = false;
            reason = "max-iterations";
            iter = 0;
            for iter = 1:obj.MaxIterations
                if norm(F, Inf) < obj.Tolerance
                    reason = "converged";
                    break;
                end
                if isempty(J) || ~useBroyden
                    J = obj.jacobian(u, F, T, z, r_cap, isCap);
                    freshJ = true;
                end

                du = obj.limitStep(obj.solveLinear(J, -F), nc, isCap);

                f0 = 0.5 * (F.' * F);
                alpha = 1.0;
                accepted = false;
                for ls = 1:obj.MaxLineSearch
                    ut = obj.project(u + alpha * du, isCap);
                    if norm(ut(1:nc)) >= obj.TrivialTol
                        [Ft, stt, okt] = obj.residual(ut, T, z, r_cap, isCap);
                        if okt && abs(stt.Z_V - stt.Z_L) >= obj.CoalescenceTol ...
                                && 0.5 * (Ft.' * Ft) <= (1 - 2 * obj.LineSearchC1 * alpha) * f0
                            accepted = true;
                            break;
                        end
                    end
                    alpha = 0.5 * alpha;
                end

                if accepted
                    if useBroyden
                        s = ut - u;
                        J = J + ((Ft - F) - J * s) * s.' / (s.' * s);
                        freshJ = false;
                    end
                    [u, F, st] = deal(ut, Ft, stt);
                elseif useBroyden && ~freshJ
                    J = [];                           % refresh with finite differences and retry
                else
                    reason = "line-search-failure";
                    break;
                end
            end
            if norm(F, Inf) < obj.Tolerance, reason = "converged"; end
            stats = obj.newtonStats(reason == "converged", reason, iter, norm(F, Inf), ...
                                    st, useBroyden, u, nc, isCap);
        end

        function [F, st, ok] = residual(obj, u, T, z, r_cap, isCap)
            nc  = numel(z);
            lnK = u(1:nc);
            PV  = exp(u(nc+1));
            if isCap
                Pcap = exp(u(nc+2));
            else
                Pcap = 0.0;
            end
            PL = obj.liquidPressure(PV, Pcap, isCap);
            xr = z .* exp(-lnK);
            S  = sum(xr);
            x  = xr / S;

            try
                [lnphiV, ZV, VV] = obj.EOS.calculateState(PV, T, z, -1, r_cap);
                [lnphiL, ZL, VL] = obj.EOS.calculateState(PL, T, x, +1, r_cap);
            catch ME
                if ~obj.isModelFailure(ME), rethrow(ME); end
                F = []; ok = false;
                st = obj.phaseState(z, x, NaN, NaN, NaN, NaN, S - 1);
                st.PL = PL;  st.Pcap = Pcap;
                return;
            end

            F = [lnK - (lnphiL + log(PL) - lnphiV - log(PV)); S - 1];
            if isCap
                Pcap_pred = obj.EOS.capillaryPressure(x, z, VL, VV, r_cap);
                F(end+1, 1) = Pcap_pred / Pcap - 1;
            end
            ok = all(isfinite(F));
            st = obj.phaseState(z, x, ZV, ZL, VV, VL, S - 1);
            st.PL = PL;  st.Pcap = Pcap;
        end

        function J = jacobian(obj, u, F0, T, z, r_cap, isCap)
            % Central differences; one-sided where a side is infeasible or
            % the model fails; zero column (handled by LM solve) if both fail.
            n = numel(u);
            J = zeros(numel(F0), n);
            for j = 1:n
                h  = obj.JacEpsilon * max(1.0, abs(u(j)));
                uf = u; uf(j) = uf(j) + h; uf = obj.project(uf, isCap);
                ub = u; ub(j) = ub(j) - h; ub = obj.project(ub, isCap);
                [Ff, ~, okf] = obj.residual(uf, T, z, r_cap, isCap);
                [Fb, ~, okb] = obj.residual(ub, T, z, r_cap, isCap);
                okf = okf && uf(j) > u(j);
                okb = okb && ub(j) < u(j);
                if okf && okb
                    J(:, j) = (Ff - Fb) / (uf(j) - ub(j));
                elseif okf
                    J(:, j) = (Ff - F0) / (uf(j) - u(j));
                elseif okb
                    J(:, j) = (F0 - Fb) / (u(j) - ub(j));
                end
            end
        end

        function du = solveLinear(obj, J, b)
            if rcond(J) > 1e-14
                du = J \ b;
            else
                JtJ = J.' * J;
                mu  = obj.JacRegularizer * max(1, norm(JtJ, 1));
                du  = (JtJ + mu * eye(size(JtJ))) \ (J.' * b);
            end
        end

        function du = limitStep(obj, du, nc, isCap)
            sK = norm(du(1:nc), Inf);
            if sK > obj.MaxStepSizeK
                du(1:nc) = du(1:nc) * (obj.MaxStepSizeK / sK);
            end
            du(nc+1) = max(min(du(nc+1), obj.MaxStepSizeP), -obj.MaxStepSizeP);
            if isCap
                du(nc+2) = max(min(du(nc+2), obj.MaxStepSizePcap), -obj.MaxStepSizePcap);
            end
        end

        function u = project(obj, u, isCap)
            nc = numel(u) - 1 - isCap;
            u(nc+1) = min(max(u(nc+1), log(obj.PressureMin)), log(obj.PressureMax));
            if isCap
                % Keep P_L = P_V - Pcap above the floor
                PV = exp(u(nc+1));
                u(nc+2) = min(u(nc+2), log(max(PV - obj.PressureFloor, 1.0)));
            end
        end

        function stats = newtonStats(obj, converged, reason, iter, resNorm, st, useBroyden, u, nc, isCap)
            if converged && norm(u(1:nc)) < obj.TrivialTol
                converged = false;
                reason = "trivial";
            end
            incipientIsLiquid = obj.identityOK(st);
            if converged && ~incipientIsLiquid
                converged = false;
                reason = "phase-identity (" + obj.PhaseIdentityTest + ")";
            end
            name = "Newton-Raphson";
            if useBroyden, name = "Quasi-Newton-Broyden"; end
            stats = struct('converged', converged, 'reason', reason, 'iterations', iter, ...
                'residual_norm', resNorm, 'solver', name + "-" + obj.modeLabel(isCap), ...
                'Z_V', st.Z_V, 'Z_L', st.Z_L, 'rho_v', 1 / st.V_V, 'rho_l', 1 / st.V_L, ...
                'sumX_res', st.sumX_res, 'rhoM_V', st.rhoM_V, 'rhoM_L', st.rhoM_L, ...
                'heavyEnriched', st.heavyEnriched, 'incipientIsLiquid', incipientIsLiquid, ...
                'Pcap', st.Pcap, 'PL', st.PL);
        end

        % ================= Helpers =================
        function st = phaseState(obj, z, x, ZV, ZL, VV, VL, sumXres)
            % Phase state at a root, with mass densities [kg/m^3] and heavy enrichment.
            MW = obj.EOS.Fluid.MW(:);
            [~, iH] = max(MW);
            st = struct('Z_V', ZV, 'Z_L', ZL, 'V_V', VV, 'V_L', VL, 'sumX_res', sumXres, ...
                'rhoM_V', 1e-3 * (z(:).' * MW) / VV, 'rhoM_L', 1e-3 * (x(:).' * MW) / VL, ...
                'heavyEnriched', x(iH) > z(iH));
        end

        function tf = identityOK(obj, st)
            massOK = st.rhoM_L > st.rhoM_V;          % false for NaN states
            compOK = st.heavyEnriched;
            switch obj.PhaseIdentityTest
                case "massDensity", tf = massOK;
                case "composition", tf = compOK;
                case "both",        tf = massOK && compOK;
                otherwise,          tf = true;
            end
        end

        function P = clampP(obj, P)
            P = min(max(P, obj.PressureMin), obj.PressureMax);
        end

        function PL = liquidPressure(obj, PV, Pcap, isCap)
            if isCap
                PL = max(PV - Pcap, obj.PressureFloor);
            else
                PL = PV;
            end
        end
    end

    methods (Static, Access = private)
        function tf = isModelFailure(ME)
            % Only EOS domain failures are recoverable; coding errors propagate.
            tf = startsWith(string(ME.identifier), "ConfinedEOS:");
        end

        function s = modeLabel(isCap)
            if isCap, s = "Combined"; else, s = "FWIOnly"; end
        end
    end
end
