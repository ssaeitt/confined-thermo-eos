classdef DewPointRunner < handle
    % DEWPOINTRUNNER Confined dew-point solving with engine caching, bulk
    % reference states, root acceptance and seeding fallbacks.
    %
    % One place for the machinery the fitting and prediction scripts share:
    %   - one ConfinedEOS / StabilityTester / FlashEngine per (mixture, rock)
    %   - bulk (r = Inf) dew point and K per case, solved once
    %   - solve() returns an ACCEPTED root only: dew closure, phase identity
    %     (mass density, from FlashEngine) and a physically possible shift
    %   - seeding routes: warm -> FWI continuation (Combined) -> bulk-seeded
    %     direct -> ConfinementScale homotopy, each with an
    %     SS -> Newton -> Broyden -> cold-Newton fallback chain.

    properties (SetAccess = private)
        WorkbookFile (1,1) string
        Engines      containers.Map        % "mixture|rock" -> struct(eos, flash, stab)
        Warm         containers.Map        % "<key>" -> struct(P, K, Pcap)
        BulkP        containers.Map        % "mixture|rock|z" -> Pa
        BulkK        containers.Map
    end

    properties (Access = public)
        PoreRadius      containers.Map     % rock name -> r_p [m]
        maxShift_psi    (1,1) double = 500
        sumXTol         (1,1) double = 1e-6
        contSteps       (1,:) double = [0.25 0.5 0.75 1]
        pGuessFactor    (1,1) double = 1.1
        combinedFromFWI (1,1) logical = true
        ssFallback      (1,1) logical = true
        FlashTolerance  (1,1) double = 1e-10
    end

    methods
        function obj = DewPointRunner(workbookFile, poreRadius, opts)
            arguments
                workbookFile (1,1) string
                poreRadius   containers.Map   % size validation is invalid here: size(Map) = [nKeys 1]
                opts.maxShift_psi (1,1) double = 500
                opts.sumXTol (1,1) double = 1e-6
                opts.contSteps (1,:) double = [0.25 0.5 0.75 1]
                opts.pGuessFactor (1,1) double = 1.1
                opts.combinedFromFWI (1,1) logical = true
                opts.ssFallback (1,1) logical = true
                opts.FlashTolerance (1,1) double = 1e-10
            end
            obj.WorkbookFile = workbookFile;
            obj.PoreRadius = poreRadius;
            fn = fieldnames(opts);
            for i = 1:numel(fn), obj.(fn{i}) = opts.(fn{i}); end
            obj.Engines = containers.Map();
            obj.Warm    = containers.Map();
            obj.BulkP   = containers.Map();
            obj.BulkK   = containers.Map();
        end

        function eng = engine(obj, mixture, rock)
            key = char(string(mixture) + "|" + string(rock));
            if ~isKey(obj.Engines, key)
                fl  = entities.FluidProperties.loadFromWorkbook(obj.WorkbookFile, string(mixture));
                rk  = entities.RockProperties.loadFromWorkbook(obj.WorkbookFile, string(rock));
                eos = thermo.ConfinedEOS(fl, rk);
                st  = solvers.StabilityTester(eos);
                fe  = solvers.FlashEngine(eos, st, 'Tolerance', obj.FlashTolerance);
                obj.Engines(key) = struct('eos', eos, 'stab', st, 'flash', fe);
            end
            eng = obj.Engines(key);
        end

        function setFWI(obj, k, p, lambda)
            % V31 ordering [k, p, lambda] -> ConfinedEOS properties, all engines.
            ev = values(obj.Engines);
            for n = 1:numel(ev)
                ev{n}.eos.k = k;  ev{n}.eos.pT_wall = p;  ev{n}.eos.lambda = lambda;
                ev{n}.eos.ConfinementScale = 1;
            end
        end

        function clearWarm(obj)
            obj.Warm = containers.Map();
        end

        function [P, K, ok] = bulk(obj, c)
            % Bulk (r = Inf) dew point and K for a case; cached. c from
            % entities.DewPointDataset.getCase.
            key = char(c.Mixture + "|" + c.Rock + "|" + strjoin(string(c.z.'), ","));
            if isKey(obj.BulkP, key)
                P = obj.BulkP(key);  K = obj.BulkK(key);  ok = isfinite(P);
                return;
            end
            eng = obj.engine(c.Mixture, c.Rock);
            kSave = eng.eos.kijc;
            eng.eos.kijc = zeros(numel(c.z));
            eng.eos.ConfinementScale = 1;
            [Pb, Kb, ~, ~, st] = eng.flash.solveDewPoint(c.T, c.PdewBulk, c.z, Inf, 'CapMode', "FWIOnly");
            eng.eos.kijc = kSave;
            ok = st.converged;
            P = NaN;  K = [];
            if ok, P = Pb;  K = Kb; end
            obj.BulkP(key) = P;  obj.BulkK(key) = K;
        end

        function res = solve(obj, c, mode, kijc, opts)
            % Confined dew point for case c in the given mode with the given
            % TIP matrix. Returns a struct with P [Pa], shift [psi], route,
            % reason, acceptance flag and root diagnostics.
            arguments
                obj
                c    (1,1) struct
                mode (1,1) string
                kijc (:,:) double
                opts.Key (1,1) string = ""          % warm-cache key; "" = auto
                opts.r   (1,1) double = NaN         % pore radius [m]; NaN = from PoreRadius
            end
            psi = entities.DewPointDataset.psi2Pa;
            eng = obj.engine(c.Mixture, c.Rock);
            eng.eos.kijc = (kijc + kijc.') / 2;
            r = opts.r;
            if ~isfinite(r), r = obj.PoreRadius(char(c.Rock)); end
            [Pb, Kb, okb] = obj.bulk(c);
            bulkRef = Pb;
            if ~okb, bulkRef = c.PdewBulk;  Kb = []; end
            key = opts.Key;
            if key == "", key = c.Mixture + "|" + c.Rock + "|" + strjoin(string(c.z.'), ",") + "|" + mode; end

            [P, Pcap, ok, info] = obj.solveConfined(eng, c, r, mode, char(key), ...
                                                    struct('P', bulkRef, 'K', Kb));
            res = info;
            res.P = P;  res.Pcap = Pcap;  res.ok = ok;  res.mode = mode;
            res.PdewBulkModel = Pb;
            res.shift_psi = (P - Pb) / psi;                 % model shift (vs model bulk)
            res.shiftExpBulk_psi = (P - c.PdewBulk) / psi;  % vs experimental bulk (V31 basis)
        end
    end

    methods (Access = private)
        function [P, Pcap, ok, info] = solveConfined(obj, eng, c, r, mode, key, bulk)
            ok = false;  P = NaN;  Pcap = NaN;
            info = solvers.DewPointRunner.emptyInfo("none");
            args = {'CapMode', mode, 'Solver', "newton"};
            eng.eos.ConfinementScale = 1;
            reasons = strings(0, 1);

            if isKey(obj.Warm, key)
                w = obj.Warm(key);
                [P, K, Pcap, st] = obj.tryFlash(eng, c, r, w.P, w.K, w.Pcap, args);
                [ok, info] = obj.acceptRoot(P, K, st, c, bulk.P);
                info.route = "warm";
                if ~ok, reasons(end+1) = "warm: " + info.reason; end
            end

            if ~ok && mode == "Combined" && obj.combinedFromFWI
                fwiKey = char(replace(string(key), "Combined", "FWIOnly"));
                if isKey(obj.Warm, fwiKey)
                    w = obj.Warm(fwiKey);
                    [P, K, Pcap, st] = obj.tryFlash(eng, c, r, w.P, w.K, [], args);
                    [ok, info] = obj.acceptRoot(P, K, st, c, bulk.P);
                    info.route = "fwi-continuation";
                    if ~ok, reasons(end+1) = "fwi-continuation: " + info.reason; end
                end
            end

            if ~ok
                [P, K, Pcap, st] = obj.tryFlash(eng, c, r, obj.pGuessFactor * bulk.P, bulk.K, [], args);
                [ok, info] = obj.acceptRoot(P, K, st, c, bulk.P);
                info.route = "direct";
                if ~ok, reasons(end+1) = "direct: " + info.reason; end
            end

            if ~ok
                [P, Pcap, ok, info] = obj.homotopy(eng, c, r, bulk, args);
                if ~ok, reasons(end+1) = info.reason; end
            end

            eng.eos.ConfinementScale = 1;
            if ok
                obj.Warm(key) = struct('P', P, 'K', exp(info.lnK), 'Pcap', Pcap);
            else
                P = NaN;  Pcap = NaN;
                info.reason = strjoin(reasons, " | ");
            end
        end

        function [P, Pcap, ok, info] = homotopy(obj, eng, c, r, bulk, args)
            P = NaN;  Pcap = NaN;  ok = false;
            info = solvers.DewPointRunner.emptyInfo("homotopy");
            steps = obj.contSteps(obj.contSteps > 0);
            Kprev = bulk.K;
            if isempty(Kprev), steps = [0, steps]; end
            Pprev = bulk.P;  Pcprev = [];
            for s = steps
                eng.eos.ConfinementScale = s;
                [Ps, Ks, Pcs, st] = obj.tryFlash(eng, c, r, Pprev, Kprev, Pcprev, args);
                [okS, infoS] = obj.acceptRoot(Ps, Ks, st, c, bulk.P);
                if ~okS
                    eng.eos.ConfinementScale = 1;
                    info.reason = "homotopy: failed at ConfinementScale = " + string(s) + ...
                                  " (" + infoS.reason + ")";
                    return;
                end
                Kprev = Ks;  Pprev = Ps;  Pcprev = Pcs;
                [P, Pcap, info] = deal(Ps, Pcs, infoS);
                info.route = "homotopy";
            end
            eng.eos.ConfinementScale = 1;
            ok = true;
        end

        function [P, K, Pcap, st] = tryFlash(obj, eng, c, r, Pguess, Kseed, Pcapseed, args)
            % Newton from the seeds; then SS -> Newton -> Broyden -> cold Newton.
            extra = {};
            if ~isempty(Kseed),    extra = [extra, {'K_seed', Kseed, 'PreconditionSS', false}]; end
            if ~isempty(Pcapseed), extra = [extra, {'Pcap_seed', Pcapseed}]; end
            [P, K, Pcap, ~, st] = eng.flash.solveDewPoint(c.T, Pguess, c.z, r, args{:}, extra{:});
            if st.converged || ~obj.ssFallback, return; end

            ssArgs = solvers.DewPointRunner.withSolver(args, "ss");
            [Pss, Kss, Pcss, ~, ~] = eng.flash.solveDewPoint(c.T, Pguess, c.z, r, ssArgs{:}, extra{:});
            if ~isfinite(Pss) || isempty(Kss) || any(~isfinite(Kss)) || norm(log(max(Kss, 1e-14))) < 0.01
                Pss = Pguess;  Kss = Kseed;  Pcss = Pcapseed;
            end
            if ~isempty(Kss)
                ex2 = {'K_seed', Kss, 'PreconditionSS', false};
                if ~isempty(Pcss) && isfinite(Pcss), ex2 = [ex2, {'Pcap_seed', Pcss}]; end
                [P2, K2, Pc2, ~, st2] = eng.flash.solveDewPoint(c.T, Pss, c.z, r, args{:}, ex2{:});
                if st2.converged || solvers.DewPointRunner.better(st2, st)
                    [P, K, Pcap, st] = deal(P2, K2, Pc2, st2);
                end
                if st.converged, return; end

                qn = solvers.DewPointRunner.withSolver(args, "quasinewton");
                [P3, K3, Pc3, ~, st3] = eng.flash.solveDewPoint(c.T, Pss, c.z, r, qn{:}, ex2{:});
                if st3.converged || solvers.DewPointRunner.better(st3, st)
                    [P, K, Pcap, st] = deal(P3, K3, Pc3, st3);
                end
                if st.converged, return; end
            end

            [P4, K4, Pc4, ~, st4] = eng.flash.solveDewPoint(c.T, Pguess, c.z, r, args{:});
            if st4.converged || solvers.DewPointRunner.better(st4, st)
                [P, K, Pcap, st] = deal(P4, K4, Pc4, st4);
            end
        end

        function [ok, info] = acceptRoot(obj, P, K, st, c, Pbulk)
            psi = entities.DewPointDataset.psi2Pa;
            info = solvers.DewPointRunner.emptyInfo("");
            if isempty(K) || any(~isfinite(K))
                info.reason = "solver: " + string(st.reason);
                ok = false;  return;
            end
            sumX = sum(c.z ./ K);
            info.sumX = sumX;  info.normLnK = norm(log(K));  info.lnK = log(K);
            info.Z_V = st.Z_V;  info.Z_L = st.Z_L;
            if isfield(st, 'rhoM_V'), info.rhoM_V = st.rhoM_V;  info.rhoM_L = st.rhoM_L; end
            info.shift_psi = (P - Pbulk) / psi;
            ok = false;
            if ~st.converged
                info.reason = "solver: " + string(st.reason);
            elseif ~isfinite(P) || P <= 0
                info.reason = "non-finite P";
            elseif abs(sumX - 1) > obj.sumXTol
                info.reason = sprintf("sum(z/K) - 1 = %.2e", sumX - 1);
            elseif ~st.incipientIsLiquid
                info.reason = sprintf("phase identity (rho_m,L = %.1f vs rho_m,V = %.1f kg/m3)", ...
                                      info.rhoM_L, info.rhoM_V);
            elseif abs(info.shift_psi) > obj.maxShift_psi
                info.reason = sprintf("|shift| = %.0f psi > %.0f psi (lost branch)", ...
                                      abs(info.shift_psi), obj.maxShift_psi);
            else
                ok = true;
                info.reason = "accepted";
            end
        end
    end

    methods (Static, Access = private)
        function info = emptyInfo(route)
            info = struct('route', route, 'reason', "", 'sumX', NaN, 'normLnK', NaN, ...
                'Z_V', NaN, 'Z_L', NaN, 'rhoM_V', NaN, 'rhoM_L', NaN, 'shift_psi', NaN, 'lnK', []);
        end

        function a = withSolver(args, solver)
            a = args;
            i = find(strcmp(a, 'Solver'), 1);
            if isempty(i), a = [a, {'Solver', solver}]; else, a{i + 1} = solver; end
        end

        function tf = better(stNew, stOld)
            tf = isfield(stNew, 'residual_norm') && isfield(stOld, 'residual_norm') && ...
                 isfinite(stNew.residual_norm) && ~(stOld.residual_norm <= stNew.residual_norm);
        end
    end
end