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
            psi = 6894.757;
            target = 100 * obj.DetectFraction;          % [% of V_sat]
            f = @(P) solvers.UncertaintyBudget.ldo(pt, P, T, z, Kdew, Vsat) - target;
            steps = [0.02, 0.1, 0.5, 2, 5, 10, 20, 40, 80, 150, 300, 600] * psi;
            lag = NaN;  prevP = Psat * (1 - 1e-7);
            for s = steps
                P = Psat - s;
                if P <= 0, break; end
                v = f(P);
                if isfinite(v) && v >= 0
                    try
                        Pd = fzero(f, [P, prevP], optimset('TolX', 0.01 * psi));
                        lag = Psat - Pd;
                    catch
                        lag = Psat - P;                     % bracket edge (conservative)
                    end
                    return;
                end
                prevP = P;
            end
        end
    end

    methods (Static)
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
