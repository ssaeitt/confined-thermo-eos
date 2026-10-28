classdef UncertaintyBudget < handle
    % UNCERTAINTYBUDGET Per-point bulk dew-point and MLDO uncertainty from the
    % laboratory accuracy figures, propagated through the bulk EOS (r_p = Inf).
    %
    % Laboratory figures (visual PVT cell, as reported by the experimental section):
    %   temperature    +- TempRel of the reading in degC (default 5 %, i.e. +-1 degC at
    %                  20 degC, +-2 at 40, +-3 at 60); dT_C overrides it when finite
    %   composition    +- CompRel relative on the heavy component(s) as charged
    %                  (default 5 % and 10 %; C1 takes up the balance)
    %   dew detection  by volume: the dew point is read when the liquid volume
    %                  reaches DetectVolume_cc in a CellVolume_cc cell
    %                  (default 0.1 cc in 20 cc, i.e. 0.5 % of V_sat)
    %   liquid volume  read to DetectVolume_cc -> MLDO reading +- 100*DetectVolume/CellVolume
    %                  points (0.5 points)
    %   pressure gauge negligible against the above (not propagated). Depletion at
    %                  1 psi/min is assumed to be at equilibrium; any lag from it acts in the
    %                  same direction as the detection lag (measured Pdew low) and is not
    %                  quantified here.
    %
    % evaluate() (dew point) returns
    %   P0     model bulk Pdew at the stated z, T          [Pa]
    %   UT     max |Pdew(T +- dT) - P0|                    [Pa]  (bound)
    %   UX     RSS over heavy components of max |Pdew(x_h (1 +- r)) - P0|,
    %          one value per entry of CompRel                [Pa]  (bound)
    %   lag    P0 - P(LDO = DetectVolume/CellVolume)         [Pa]  (one-sided)
    %          the measured value lies BELOW the true dew point by up to lag
    %
    % evaluateMLDO() returns the same budget for the CCE maximum liquid dropout
    %   M0, UT, UX(q) [points] and Uread (reading resolution, points).
    %
    % Bands and weights (static):
    %   band  : model Pdew acceptable in [Pm - U, Pm + lag + U], U = sqrt(UX^2 + UT^2 + Urep^2)
    %   sigma : rectangular bounds -> standard uncertainties
    %           sigma^2 = (UX^2 + UT^2 + Urep^2) / 3 + lag^2 / 12, target Pm + lag / 2
    %
    % GUM / NIST TN 1297 combination (gum, gumMLDO):
    %   Each excursion above is a Type B input. With the distribution assigned to
    %   it (CompDist, TempDist, GaugeDist, ReadDist) it becomes a standard
    %   uncertainty u = bound / divisor:
    %       "rectangular" sqrt(3) | "triangular" sqrt(6) | "normal95" 2 (bound given
    %       at ~95 %) | "standard" 1 (the figure already is a standard uncertainty)
    %   One-sided effects (detection lag, dew pick during the pressure ramp,
    %   width w) lie on [0, w]. The dew pick is rectangular: correction +w/2,
    %   u = w/sqrt(12). The detection lag follows LagDist:
    %       "rectangular" correction +w/2, u = w/sqrt(12) (liquid read at the
    %                     volume DetectVolume_cc, anywhere up to w)
    %       "cloudpoint"  right-triangular with its mode at 0: correction +w/3,
    %                     u = w/sqrt(18). The dew point is taken at the first
    %                     cloud-like appearance (Doryani thesis, Sec. 3.2), so a
    %                     lag near zero is the most likely value and w (the liquid
    %                     volume DetectVolume_cc) is its upper bound.
    %   Combined standard uncertainty (law of propagation, uncorrelated inputs,
    %   sensitivities by finite difference):
    %       u_c^2 = u_X^2 + u_T^2 + u_gauge^2 + u_lag^2 + u_pick^2 + u_rep^2
    %   Expanded uncertainty U_c = k u_c (Coverage k, default 2: ~0.95 level of
    %   confidence). Model check: E_n = (model - y) / U_c, |E_n| <= 1 consistent.
    %   u_rep is a Type A (repeatability) standard uncertainty passed by the caller.
    %
    % The sensitivities are evaluated with the BIP matrix of the fluid passed in.

    properties (SetAccess = private)
        CompRel          (1,:) double = [0.05, 0.10]
        TempRel          (1,1) double = 0.05
        dT_C             (1,1) double = NaN
        CellVolume_cc    (1,1) double = 20
        DetectVolume_cc  (1,1) double = 0.1
        LightComponent   (1,1) string = "C1"
        PminFrac         (1,1) double = 0.10       % CCE for MLDO down to PminFrac * Psat
        nCoarse          (1,1) double = 60         % CCE points below 0.95 Psat
        % GUM settings
        CompDist         (1,1) string = "rectangular"
        TempDist         (1,1) string = "rectangular"
        ReadDist         (1,1) string = "rectangular"
        GaugeDist        (1,1) string = "rectangular"
        LagDist          (1,1) string = "rectangular"   % "rectangular" | "cloudpoint"
        Gauge_psi        (1,1) double = 0          % pressure transducer offset bound (or u if "standard")
        GaugeRel         (1,1) double = 0          % pressure transducer span bound, fraction of reading
        PickWidth_psi    (1,1) double = 0          % one-sided: ramp rate x observation interval
        Coverage         (1,1) double = 2          % k
    end

    properties (Dependent)
        DetectFraction   % DetectVolume_cc / CellVolume_cc
        Uread            % MLDO reading resolution [points]
    end

    methods
        function obj = UncertaintyBudget(opts)
            arguments
                opts.CompRel (1,:) double = [0.05, 0.10]
                opts.TempRel (1,1) double = 0.05
                opts.dT_C (1,1) double = NaN
                opts.CellVolume_cc (1,1) double = 20
                opts.DetectVolume_cc (1,1) double = 0.1
                opts.LightComponent (1,1) string = "C1"
                opts.PminFrac (1,1) double = 0.10
                opts.nCoarse (1,1) double = 60
                opts.CompDist (1,1) string {mustBeMember(opts.CompDist, ["rectangular", "triangular", "normal95", "standard"])} = "rectangular"
                opts.TempDist (1,1) string {mustBeMember(opts.TempDist, ["rectangular", "triangular", "normal95", "standard"])} = "rectangular"
                opts.ReadDist (1,1) string {mustBeMember(opts.ReadDist, ["rectangular", "triangular", "normal95", "standard"])} = "rectangular"
                opts.GaugeDist (1,1) string {mustBeMember(opts.GaugeDist, ["rectangular", "triangular", "normal95", "standard"])} = "rectangular"
                opts.LagDist (1,1) string {mustBeMember(opts.LagDist, ["rectangular", "cloudpoint"])} = "rectangular"
                opts.Gauge_psi (1,1) double = 0
                opts.GaugeRel (1,1) double = 0
                opts.PickWidth_psi (1,1) double = 0
                opts.Coverage (1,1) double = 2
            end
            fn = fieldnames(opts);
            for i = 1:numel(fn), obj.(fn{i}) = opts.(fn{i}); end
        end

        function f = get.DetectFraction(obj), f = obj.DetectVolume_cc / obj.CellVolume_cc; end
        function u = get.Uread(obj),          u = 100 * obj.DetectFraction; end

        function dT = tempUncertainty(obj, T)
            % Temperature bound [K] at T [K]: dT_C if set, else TempRel of the reading in degC.
            if isfinite(obj.dT_C), dT = obj.dT_C; else, dT = obj.TempRel * (T - 273.15); end
        end

        function u = evaluate(obj, fluid, rock, z, T, Pmeas)
            % fluid : entities.FluidProperties (BIP matrix as wanted)
            % rock  : entities.RockProperties (only to build the EOS; no effect at r = Inf)
            % z, T [K], Pmeas [Pa]
            z = z(:) / sum(z);
            eos = thermo.ConfinedEOS(fluid, rock);
            sat = solvers.SaturationSolver(eos);
            pt  = solvers.PTFlash(eos);
            u = struct('P0', NaN, 'dT', obj.tempUncertainty(T), 'UT', NaN, 'dP_T', [NaN, NaN], ...
                       'UX', nan(1, numel(obj.CompRel)), 'UX_parts', {cell(1, numel(obj.CompRel))}, ...
                       'oneSided', false(1, numel(obj.CompRel)), 'lag', NaN, 'Vsat', NaN, 'note', "");

            [r0, br] = solvers.UncertaintyBudget.upperDew(sat, pt, T, Pmeas, z);
            if ~r0.converged
                u.note = "bulk dew point: " + r0.reason;  return;
            end
            u.P0 = r0.P;
            if br ~= "upper dew", u.note = br; end

            % Temperature
            for s = [-1, 1]
                rT = solvers.UncertaintyBudget.upperDew(sat, pt, T + s * u.dT, u.P0, z);
                if rT.converged, u.dP_T((s + 3) / 2) = rT.P - u.P0; end
            end
            u.UT = max(abs(u.dP_T), [], 'omitnan');

            % Composition: each heavy component separately, C1 takes up the balance
            iL = fluid.componentIndex(obj.LightComponent);
            heavy = setdiff(1:fluid.NC, iL);
            for q = 1:numel(obj.CompRel)
                parts = nan(1, numel(heavy));
                for h = 1:numel(heavy)
                    d = [NaN, NaN];
                    for s = [-1, 1]
                        zz = obj.perturb(z, iL, heavy(h), s * obj.CompRel(q));
                        if isempty(zz), continue; end
                        rX = solvers.UncertaintyBudget.upperDew(sat, pt, T, u.P0, zz);
                        if rX.converged, d((s + 3) / 2) = rX.P - u.P0; end
                    end
                    parts(h) = max(abs(d), [], 'omitnan');
                    if any(isnan(d)), u.oneSided(q) = true; end
                end
                u.UX_parts{q} = parts;
                u.UX(q) = sqrt(sum(parts.^2, 'omitnan'));
            end

            % Detection lag: liquid volume reaches DetectFraction of V_sat
            [~, ~, u.Vsat] = eos.calculateState(u.P0, T, z, -1, Inf);
            u.lag = obj.detectionLag(pt, T, z, u.P0, r0.K, u.Vsat);
        end

        function m = evaluateMLDO(obj, fluid, rock, z, T, Pmeas)
            % CCE maximum liquid dropout [% of V_sat] and its budget [points].
            % Each perturbed state gets its own upper dew point and V_sat (as in the
            % laboratory, where V_sat is read for the fluid actually charged).
            z = z(:) / sum(z);
            eos = thermo.ConfinedEOS(fluid, rock);
            sat = solvers.SaturationSolver(eos);
            pt  = solvers.PTFlash(eos);
            pvt = solvers.PVTExperiments(pt);
            nQ  = numel(obj.CompRel);
            m = struct('M0', NaN, 'P0', NaN, 'dT', obj.tempUncertainty(T), 'dM_T', [NaN, NaN], 'UT', NaN, ...
                       'UX', nan(1, nQ), 'UX_parts', {cell(1, nQ)}, 'Uread', obj.Uread, 'note', "");
            [m.M0, m.P0, m.note] = obj.mldoAt(sat, pt, pvt, T, Pmeas, z);
            if ~isfinite(m.M0), return; end

            for s = [-1, 1]
                m.dM_T((s + 3) / 2) = obj.mldoAt(sat, pt, pvt, T + s * m.dT, m.P0, z) - m.M0;
            end
            m.UT = max(abs(m.dM_T), [], 'omitnan');

            iL = fluid.componentIndex(obj.LightComponent);
            heavy = setdiff(1:fluid.NC, iL);
            for q = 1:nQ
                parts = nan(1, numel(heavy));
                for h = 1:numel(heavy)
                    d = [NaN, NaN];
                    for s = [-1, 1]
                        zz = obj.perturb(z, iL, heavy(h), s * obj.CompRel(q));
                        if isempty(zz), continue; end
                        d((s + 3) / 2) = obj.mldoAt(sat, pt, pvt, T, m.P0, zz) - m.M0;
                    end
                    parts(h) = max(abs(d), [], 'omitnan');
                end
                m.UX_parts{q} = parts;
                m.UX(q) = sqrt(sum(parts.^2, 'omitnan'));
            end
        end
    end

    methods
        function g = gum(obj, P0, Pm, UX, UT, lag, uRep)
            % GUM combination for one dew point. All pressures in the same unit
            % (Pa as returned by evaluate). P0 model, Pm measured, UX (1 x nQ) and
            % UT the excursions of evaluate, lag the detection lag, uRep Type A.
            if nargin < 7 || ~isfinite(uRep), uRep = 0; end
            psi  = 6894.757;
            if ~isfinite(lag), lag = 0; end
            pick = obj.PickWidth_psi * psi;
            [cLag, dLag] = obj.lagFactors();
            g.y     = Pm + cLag * lag + pick / 2;            % best estimate of the measurand
            g.uX    = UX / obj.divisor(obj.CompDist);
            g.uT    = UT / obj.divisor(obj.TempDist);
            g.uGoff  = obj.Gauge_psi * psi / obj.divisor(obj.GaugeDist);   % offset part
            g.uGspan = obj.GaugeRel * Pm / obj.divisor(obj.GaugeDist);     % span part (at this reading)
            g.uG     = hypot(g.uGoff, g.uGspan);
            g.uLag  = lag / dLag;
            g.uPick = pick / sqrt(12);
            g.uRep  = uRep;
            g.uTypeB = sqrt(g.uX.^2 + g.uT^2 + g.uG^2 + g.uLag^2 + g.uPick^2);
            g.uc    = sqrt(g.uTypeB.^2 + g.uRep^2);
            g.Uc    = obj.Coverage * g.uc;
            g.dev   = P0 - g.y;
            g.En    = g.dev ./ g.Uc;
        end

        function g = gumMLDO(obj, M0, Mm, UX, UT)
            % GUM combination for the CCE maximum liquid dropout [points].
            g.y     = Mm;
            g.uX    = UX / obj.divisor(obj.CompDist);
            g.uT    = UT / obj.divisor(obj.TempDist);
            g.uRead = obj.Uread / obj.divisor(obj.ReadDist);
            g.uc    = sqrt(g.uX.^2 + g.uT^2 + g.uRead^2);
            g.Uc    = obj.Coverage * g.uc;
            g.dev   = M0 - g.y;
            g.En    = g.dev ./ g.Uc;
        end

        function [c, d] = lagFactors(obj)
            % Detection lag of width w: correction c*w, standard uncertainty w/d.
            switch obj.LagDist
                case "cloudpoint", c = 1/3;  d = sqrt(18);
                otherwise,         c = 1/2;  d = sqrt(12);
            end
        end

        function u = uTemp(obj, T)
            % Standard uncertainty of the temperature [K] at T [K].
            u = obj.tempUncertainty(T) / obj.divisor(obj.TempDist);
        end

        function u = uCompRel(obj)
            % Relative standard uncertainty of the heavy fraction(s), one per CompRel.
            u = obj.CompRel / obj.divisor(obj.CompDist);
        end
    end

    methods (Access = private)
        function [M, P, note] = mldoAt(obj, sat, pt, pvt, T, Pguess, z)
            [M, P, note] = deal(NaN, NaN, "");
            [r, br] = solvers.UncertaintyBudget.upperDew(sat, pt, T, Pguess, z);
            if ~r.converged, note = "dew: " + r.reason; return; end
            if br ~= "upper dew", note = br; end
            P = r.P;
            tbl = pvt.cce(T, z, Inf, P, solvers.PVTExperiments.pressureGrid(P, obj.PminFrac * P, obj.nCoarse));
            M = solvers.PVTExperiments.mldo(tbl);
        end

        function lag = detectionLag(obj, pt, T, z, Psat, Kdew, Vsat)
            % Pressure drop below Psat at which the liquid volume reaches
            % DetectFraction of V_sat. Coarse scan for the first step with enough
            % liquid, then bisection between the last converged point without
            % enough liquid (Psat itself to start with: no liquid at the dew point)
            % and that step. Flash failures are skipped, never used as bracket ends
            % (the former fzero bracket could take a failed point and fell back to
            % the step edge, so the lag came out as one of the step sizes).
            psi = 6894.757;
            target = 100 * obj.DetectFraction;          % [% of V_sat]
            f = @(P) solvers.UncertaintyBudget.ldo(pt, P, T, z, Kdew, Vsat) - target;
            steps = [0.02, 0.1, 0.5, 1, 2, 5, 10, 20, 40, 80, 150, 300, 600] * psi;
            lag = NaN;  Plo = Psat;
            for s = steps
                P = Psat - s;
                if P <= 0, break; end
                v = f(P);
                if ~isfinite(v), continue; end
                if v < 0, Plo = P; continue; end
                Phi = P;
                for it = 1:60
                    if Plo - Phi < 0.01 * psi, break; end
                    Pm = 0.5 * (Plo + Phi);  vm = f(Pm);
                    if ~isfinite(vm)                     % flash failed: try off-centre
                        Pm = Plo - 0.3 * (Plo - Phi);  vm = f(Pm);
                        if ~isfinite(vm), break; end
                    end
                    if vm >= 0, Phi = Pm; else, Plo = Pm; end
                end
                lag = Psat - 0.5 * (Plo + Phi);
                return;
            end
        end
    end

    methods (Static)
        function d = divisor(dist)
            % Bound -> standard uncertainty (GUM 4.3.7, 4.3.9; NIST TN 1297 4.6).
            switch string(dist)
                case "rectangular", d = sqrt(3);
                case "triangular",  d = sqrt(6);
                case "normal95",    d = 2;
                case "standard",    d = 1;
                otherwise, error('UncertaintyBudget:Distribution', 'Unknown distribution "%s".', dist);
            end
        end

        function s = report(y, U)
            % "y +- U" with U to two significant digits and y to the same decimal
            % place (no more than one digit beyond the uncertainty; JCED rule).
            if ~isfinite(U) || U <= 0 || ~isfinite(y), s = string(sprintf('%g', y)); return; end
            d  = 1 - floor(log10(U));
            Ur = round(U, d);
            if Ur >= 10^(2 - d), d = d - 1; Ur = round(U, d); end     % 9.96 -> 10
            yr = round(y, d);
            fmt = sprintf('%%.%df', max(d, 0));
            s = string(sprintf([fmt ' ' char(177) ' ' fmt], yr, Ur));
        end

        function zz = perturb(z, iL, ih, rel)
            % Heavy component ih scaled by (1 + rel); the light component takes up
            % the balance. Empty when the light fraction would become non-positive.
            zz = z;  zz(ih) = z(ih) * (1 + rel);
            zz(iL) = 1 - (sum(zz) - zz(iL));
            if zz(iL) <= 0, zz = []; end
        end

        function v = ldo(pt, P, T, z, Kseed, Vsat)
            % Liquid volume in % of V_sat at P (CCE basis, PVTExperiments.cce).
            r = pt.solve(P, T, z, Inf, 'K_seed', Kseed);
            if ~r.converged, r = pt.solve(P, T, z, Inf); end     % retry without the warm start
            if ~r.converged, v = NaN; return; end
            if r.phases == 1, v = 0; return; end
            v = 100 * r.SL * r.V_total / Vsat;
        end

        function [lo, hi, U] = band(Pm, UX, UT, lag, Urep)
            % Acceptance band for the model Pdew [same units as inputs]. Urep: optional
            % repeatability bound added in quadrature (default 0).
            if nargin < 5, Urep = 0; end
            U  = sqrt(UX.^2 + UT.^2 + Urep.^2);
            lo = Pm - U;
            hi = Pm + lag + U;
        end

        function [target, sigma] = weight(Pm, UX, UT, lag, Urep)
            % Target and standard uncertainty for a weighted least-squares fit.
            if nargin < 5, Urep = 0; end
            target = Pm + lag / 2;
            sigma  = sqrt((UX.^2 + UT.^2 + Urep.^2) / 3 + lag.^2 / 12);
        end

        function [res, branch] = upperDew(sat, pt, T, Pguess, z)
            % Dew point on the UPPER (retrograde) branch: single phase just above,
            % two phases just below (as MLDO_check localUpperDew).
            branch = "not checked";
            res = struct('converged', false, 'reason', "not run", 'P', NaN, 'K', []);
            for g = [1, 0.97, 1.03, 1.3, 0.8, 1.6]
                res = sat.solve("dew", T, g * Pguess, z, Inf);
                if ~res.converged, continue; end
                above = pt.solve(res.P * 1.002, T, z, Inf);
                below = pt.solve(res.P * 0.998, T, z, Inf, 'K_seed', res.K);
                if above.phases == 1 && below.phases == 2
                    branch = "upper dew";
                    return;
                end
                branch = "lower dew branch rejected";
            end
            if res.converged && branch ~= "upper dew"
                res.converged = false;
                res.reason = "only a lower-branch dew point was found";
            end
        end
    end
end
