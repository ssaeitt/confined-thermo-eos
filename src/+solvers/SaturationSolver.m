classdef SaturationSolver < handle
    % SATURATIONSOLVER Bubble- and dew-point pressure of confined mixtures,
    % FWI-only, in one incipient-phase formulation.
    %
    %   Unknowns  u = [ln K (nc); ln P],   K_i = y_i / x_i
    %   Dew   : feed = vapour (y = z),  incipient liquid  x = (z ./ K) / sum(z ./ K)
    %           F_i = ln K_i - [ln phi_i^L(x) - ln phi_i^V(z)],  F_nc+1 = sum(z ./ K) - 1
    %   Bubble: feed = liquid (x = z),  incipient vapour  y = (z .* K) / sum(z .* K)
    %           F_i = ln K_i - [ln phi_i^L(z) - ln phi_i^V(y)],  F_nc+1 = sum(z .* K) - 1
    % (ln P cancels: one pressure for both phases in FWI-only.)
    %
    % Roots are forced (+1 liquid, -1 vapour), as in FlashEngine, to keep the
    % two phases distinct; trivial / coalesced states are rejected and the
    % phase identity is checked by MASS density.
    %
    % Seeding, in order: user K_seed -> two-sided TPD at P_guess and at a
    % scan of lower pressures (PhaseStability) -> Wilson. If the confined solve fails, a ConfinementScale homotopy
    % from the bulk solution is attempted (bulk -> 0.25 -> ... -> 1), the
    % same strategy DewPointRunner uses for dew points.
    %
    % Dew points of the frozen fitting pipeline should still come from
    % DewPointRunner/FlashEngine; this class is the general-purpose solver
    % for envelopes, bubble points and PVT experiments, and is verified
    % against FlashEngine on the dew side (see PVT_verification.m).
    %
    % Note: a lean feed has no bubble point at a given T (it lies right of
    % the critical point); a failed bubble solve on such a feed is physics,
    % not a solver defect. Dew curves can also have an upper and a lower
    % branch; the branch found depends on P_guess.

    properties (SetAccess = private)
        EOS       thermo.ConfinedEOS
        Stability solvers.PhaseStability

        MaxIterations  (1,1) double {mustBeInteger, mustBePositive} = 100
        Tolerance      (1,1) double {mustBePositive} = 1e-10
        JacEpsilon     (1,1) double {mustBePositive} = 1e-6
        MaxStepLnK     (1,1) double {mustBePositive} = 0.20
        MaxStepLnP     (1,1) double {mustBePositive} = 0.05
        LineSearchC1   (1,1) double {mustBePositive} = 1e-4
        MaxLineSearch  (1,1) double {mustBeInteger, mustBePositive} = 14
        TrivialTol     (1,1) double {mustBePositive} = 0.01
        CoalescenceTol (1,1) double {mustBePositive} = 1e-5
        PressureMin    (1,1) double {mustBePositive} = 1e5
        PressureMax    (1,1) double {mustBePositive} = 35000 * 6894.757
        ContSteps      (1,:) double = [0.25 0.5 0.75 1]
        ScanFactors    (1,:) double {mustBePositive} = [1 0.97 0.93 0.88 0.8 0.7 0.55]
        NonmonotoneMemory (1,1) double {mustBeInteger, mustBePositive} = 5
        MeritFloor     (1,1) double {mustBePositive} = 1e-8   % 0.5||F||^2 below which growth is tolerated
    end

    methods
        function obj = SaturationSolver(eosEngine, opts)
            arguments
                eosEngine (1,1) thermo.ConfinedEOS
                opts.Stability solvers.PhaseStability = solvers.PhaseStability.empty
                opts.MaxIterations  (1,1) double
                opts.Tolerance      (1,1) double
                opts.JacEpsilon     (1,1) double
                opts.MaxStepLnK     (1,1) double
                opts.MaxStepLnP     (1,1) double
                opts.LineSearchC1   (1,1) double
                opts.MaxLineSearch  (1,1) double
                opts.TrivialTol     (1,1) double
                opts.CoalescenceTol (1,1) double
                opts.PressureMin    (1,1) double
                opts.PressureMax    (1,1) double
                opts.ContSteps      (1,:) double
                opts.ScanFactors    (1,:) double
                opts.NonmonotoneMemory (1,1) double
                opts.MeritFloor     (1,1) double
            end
            obj.EOS = eosEngine;
            if isempty(opts.Stability)
                obj.Stability = solvers.PhaseStability(eosEngine);
            else
                if opts.Stability.EOS ~= eosEngine
                    error('SaturationSolver:EngineMismatch', ...
                          'PhaseStability must wrap the same ConfinedEOS instance.');
                end
                obj.Stability = opts.Stability;
            end
            opts = rmfield(opts, 'Stability');
            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
        end

        function res = solve(obj, type, T, P_guess, z, r_cap, opts)
            % type "dew" | "bubble"; T [K]; P_guess [Pa]; z feed; r_cap [m].
            % opts.K_seed      : K = y/x seed (e.g. from a nearby flash or the
            %                    bulk solution).
            % opts.Continuation: allow the bulk -> confined homotopy fallback.
            arguments
                obj
                type    (1,1) string {mustBeMember(type, ["dew", "bubble"])}
                T       (1,1) double {mustBePositive}
                P_guess (1,1) double {mustBePositive}
                z       (:,1) double {mustBeNonnegative}
                r_cap   (1,1) double {mustBePositive}
                opts.K_seed double = []
                opts.Continuation (1,1) logical = true
            end
            nc = obj.EOS.Fluid.NC;
            if numel(z) ~= nc
                error('SaturationSolver:CompositionSize', 'z has %d entries; NC = %d.', numel(z), nc);
            end
            z = max(z, 1e-20);
            z = z / sum(z);

            res = obj.direct(type, T, P_guess, z, r_cap, opts.K_seed(:));
            res.route = "direct";
            if res.converged || ~opts.Continuation || ~isfinite(r_cap)
                return;
            end

            resH = obj.homotopy(type, T, P_guess, z, r_cap);
            if resH.converged
                res = resH;
            else
                res.reason = "direct: " + res.reason + " | " + resH.reason;
            end
        end
    end

    methods (Access = private)
        function res = direct(obj, type, T, P_guess, z, r_cap, Kseed)
            % Newton from each seed candidate in turn; first accepted root wins.
            %   1. user K_seed at P_guess
            %   2. two-sided TPD at P_guess * ScanFactors (stop at the first
            %      unstable state whose trial phase is the incipient type)
            %   3. Wilson at P_guess
            % Near-critical feeds make the TPD stationary point at P_guess
            % almost trivial; scanning into the two-phase region gives a
            % well-separated seed (the failure mode that motivated this).
            res = [];
            if ~isempty(Kseed)
                res = obj.newton(type, T, P_guess, z, r_cap, Kseed);
                res.seed = "user";
                if res.converged, return; end
            end
            want = "liquid";
            if type == "bubble", want = "vapor"; end
            for f = obj.ScanFactors
                Pf = obj.clampP(f * P_guess);
                try
                    st = obj.Stability.test(Pf, T, z, r_cap);
                catch ME
                    if ~startsWith(string(ME.identifier), "ConfinedEOS:"), rethrow(ME); end
                    continue;
                end
                if st.isUnstable && st.trialType == want && norm(log(st.K)) > 0.05
                    r = obj.newton(type, T, Pf, z, r_cap, st.K);
                    r.seed = "tpd@" + string(f);
                    if r.converged, res = r; return; end
                    if isempty(res), res = r; end
                end
            end
            r = obj.newton(type, T, P_guess, z, r_cap, obj.EOS.wilsonK(P_guess, T));
            r.seed = "wilson";
            if r.converged || isempty(res)
                res = r;
            end
        end

        function res = homotopy(obj, type, T, P_guess, z, r_cap)
            scale0 = obj.EOS.ConfinementScale;
            restore = onCleanup(@() obj.resetScale(scale0)); %#ok<NASGU>

            obj.EOS.ConfinementScale = 1;
            rb = obj.direct(type, T, P_guess, z, Inf, []);
            if ~rb.converged
                res = rb;
                res.reason = "homotopy: bulk solve failed (" + rb.reason + ")";
                res.route = "homotopy";
                return;
            end
            P = rb.P;  K = rb.K;
            steps = obj.ContSteps(obj.ContSteps > 0 & obj.ContSteps <= 1);
            if isempty(steps) || steps(end) ~= 1
                steps = [steps, 1];                    % always finish at full confinement
            end
            for s = steps
                obj.EOS.ConfinementScale = s;
                rs = obj.direct(type, T, P, z, r_cap, K);
                if ~rs.converged
                    res = rs;
                    res.reason = "homotopy: failed at ConfinementScale = " + string(s) + ...
                                 " (" + rs.reason + ")";
                    res.route = "homotopy";
                    return;
                end
                P = rs.P;  K = rs.K;
            end
            res = rs;
            res.route = "homotopy";
            res.seed = "bulk";
        end

        function resetScale(obj, s)
            obj.EOS.ConfinementScale = s;
        end

        function res = newton(obj, type, T, P0, z, r_cap, K0)
            nc = numel(z);
            u = [log(max(K0(:), 1e-14)); log(obj.clampP(P0))];
            [F, st, ok] = obj.residual(type, u, T, z, r_cap);
            reason = "max-iterations";
            it = 0;
            if ~ok
                res = obj.pack(type, u, st, false, "initial-state-failure", 0, z);
                return;
            end
            meritHist = 0.5 * (F.' * F);
            for it = 1:obj.MaxIterations
                if norm(F, Inf) < obj.Tolerance
                    reason = "converged";
                    break;
                end
                J = zeros(nc + 1);
                for j = 1:nc + 1
                    h  = obj.JacEpsilon * max(1, abs(u(j)));
                    up = u; up(j) = up(j) + h;
                    um = u; um(j) = um(j) - h;
                    [Fp, ~, okp] = obj.residual(type, up, T, z, r_cap);
                    [Fm, ~, okm] = obj.residual(type, um, T, z, r_cap);
                    if okp && okm
                        J(:, j) = (Fp - Fm) / (2 * h);
                    elseif okp
                        J(:, j) = (Fp - F) / h;
                    elseif okm
                        J(:, j) = (F - Fm) / h;
                    end
                end
                if rcond(J) > 1e-14
                    du = -(J \ F);
                else
                    du = -((J.' * J + 1e-10 * max(1, norm(J, 1)^2) * eye(nc + 1)) \ (J.' * F));
                end
                sK = norm(du(1:nc), Inf);
                if sK > obj.MaxStepLnK
                    du(1:nc) = du(1:nc) * (obj.MaxStepLnK / sK);
                end
                du(nc+1) = max(min(du(nc+1), obj.MaxStepLnP), -obj.MaxStepLnP);

                % Non-monotone (Grippo-Lampariello-Lucidi) reference with an
                % absolute floor: near the critical region the full Newton
                % step can raise ||F|| transiently while converging, and a
                % monotone Armijo test then crawls at alpha ~ 1/64.
                fRef = max(max(meritHist(max(1, end - obj.NonmonotoneMemory + 1):end)), obj.MeritFloor);
                alpha = 1.0;
                accepted = false;
                for ls = 1:obj.MaxLineSearch
                    ut = u + alpha * du;
                    ut(nc+1) = log(obj.clampP(exp(ut(nc+1))));
                    if norm(ut(1:nc)) >= obj.TrivialTol
                        [Ft, stt, okt] = obj.residual(type, ut, T, z, r_cap);
                        if okt && abs(stt.Z_V - stt.Z_L) >= obj.CoalescenceTol ...
                                && 0.5 * (Ft.' * Ft) <= (1 - 2 * obj.LineSearchC1 * alpha) * fRef
                            accepted = true;
                            break;
                        end
                    end
                    alpha = 0.5 * alpha;
                end
                if ~accepted
                    reason = "line-search-failure";
                    break;
                end
                [u, F, st] = deal(ut, Ft, stt);
                meritHist(end+1) = 0.5 * (F.' * F); %#ok<AGROW>
            end
            if norm(F, Inf) < obj.Tolerance, reason = "converged"; end

            converged = reason == "converged";
            if converged && norm(u(1:nc)) < obj.TrivialTol
                converged = false;  reason = "trivial";
            elseif converged && ~st.identityOK
                converged = false;
                reason = sprintf("phase identity (rho_m,L = %.1f, rho_m,V = %.1f kg/m3)", st.rhoM_L, st.rhoM_V);
            end
            res = obj.pack(type, u, st, converged, reason, it, z);
            res.residual_norm = norm(F, Inf);
        end

        function [F, st, ok] = residual(obj, type, u, T, z, r_cap)
            nc  = numel(z);
            lnK = u(1:nc);
            P   = exp(u(nc+1));
            MW  = obj.EOS.Fluid.MW(:);
            st  = struct('Z_V', NaN, 'Z_L', NaN, 'V_V', NaN, 'V_L', NaN, 'x', z, 'y', z, ...
                         'rhoM_L', NaN, 'rhoM_V', NaN, 'identityOK', false);
            if type == "dew"
                wr = z .* exp(-lnK);  S = sum(wr);  x = wr / S;  y = z;
            else
                wr = z .* exp(lnK);   S = sum(wr);  y = wr / S;  x = z;
            end
            try
                [lnphiL, ZL, VL] = obj.EOS.calculateState(P, T, x, +1, r_cap);
                [lnphiV, ZV, VV] = obj.EOS.calculateState(P, T, y, -1, r_cap);
            catch ME
                if ~startsWith(string(ME.identifier), "ConfinedEOS:"), rethrow(ME); end
                F = [];  ok = false;
                return;
            end
            F  = [lnK - (lnphiL - lnphiV); S - 1];
            ok = all(isfinite(F));
            rhoL = 1e-3 * (x.' * MW) / VL;
            rhoV = 1e-3 * (y.' * MW) / VV;
            st = struct('Z_V', ZV, 'Z_L', ZL, 'V_V', VV, 'V_L', VL, 'x', x, 'y', y, ...
                        'rhoM_L', rhoL, 'rhoM_V', rhoV, 'identityOK', rhoL > rhoV);
        end

        function res = pack(~, type, u, st, converged, reason, it, z)
            nc = numel(z);
            if type == "dew"
                incipient = st.x;
            else
                incipient = st.y;
            end
            res = struct('type', type, 'P', exp(u(nc+1)), 'K', exp(u(1:nc)), ...
                'converged', converged, 'reason', string(reason), 'iterations', it, ...
                'residual_norm', NaN, 'incipient', incipient, 'x', st.x, 'y', st.y, ...
                'Z_L', st.Z_L, 'Z_V', st.Z_V, 'V_L', st.V_L, 'V_V', st.V_V, ...
                'rhoM_L', st.rhoM_L, 'rhoM_V', st.rhoM_V, 'seed', "", 'route', "");
        end

        function P = clampP(obj, P)
            P = min(max(P, obj.PressureMin), obj.PressureMax);
        end
    end
end
