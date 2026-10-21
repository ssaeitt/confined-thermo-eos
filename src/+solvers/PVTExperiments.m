classdef PVTExperiments < handle
    % PVTEXPERIMENTS Constant-composition expansion (CCE) and constant-volume
    % depletion (CVD) built on PTFlash, FWI-only, single pore radius.
    %
    % Volumes are Peneloux-translated molar volumes from ConfinedEOS.
    % Liquid dropout is reported as the liquid volume fraction:
    %   CCE : SL = V_liq / V_total(P)
    %   CVD : SL = V_liq / V_cell     (V_cell = feed volume at P_sat)
    %
    % Caveat for the confined model: the EOS is calibrated on Pdew only.
    % Dropout, densities and produced-gas compositions below Psat are model
    % predictions, best reported as confined-vs-bulk differences.

    properties (SetAccess = private)
        EOS   thermo.ConfinedEOS
        Flash solvers.PTFlash
    end

    methods
        function obj = PVTExperiments(flashEngine)
            arguments
                flashEngine (1,1) solvers.PTFlash
            end
            obj.Flash = flashEngine;
            obj.EOS   = flashEngine.EOS;
        end

        function tbl = cce(obj, T, z, r_cap, Psat, Plist, opts)
            % Psat [Pa]: saturation pressure (from SaturationSolver or
            % DewPointRunner); Plist [Pa]: pressures, any order.
            % opts.SatType "dew" | "bubble" sets which root defines V_sat.
            arguments
                obj
                T     (1,1) double {mustBePositive}
                z     (:,1) double {mustBeNonnegative}
                r_cap (1,1) double {mustBePositive}
                Psat  (1,1) double {mustBePositive}
                Plist (:,1) double {mustBePositive}
                opts.SatType (1,1) string {mustBeMember(opts.SatType, ["dew", "bubble"])} = "dew"
            end
            z = z / sum(z);
            flag = -1;
            if opts.SatType == "bubble", flag = +1; end
            [~, ~, Vsat] = obj.EOS.calculateState(Psat, T, z, flag, r_cap);

            Plist = sort(Plist, 'descend');
            n = numel(Plist);
            [Vrel, SL, beta, rhoL, rhoV] = deal(nan(n, 1));
            phases = zeros(n, 1);
            reason = strings(n, 1);
            K = [];
            for i = 1:n
                r = obj.Flash.solve(Plist(i), T, z, r_cap, 'K_seed', K);
                phases(i) = r.phases;  reason(i) = r.reason;
                if ~r.converged
                    continue;                                   % row stays NaN, reason kept
                end
                Vrel(i) = r.V_total / Vsat;
                SL(i)   = r.SL;
                beta(i) = r.beta;
                rhoL(i) = r.rhoM_L;  rhoV(i) = r.rhoM_V;
                if r.phases == 2, K = r.K; end
            end
            tbl = table(Plist, Plist / 6894.757, phases, beta, Vrel, SL, rhoL, rhoV, reason, ...
                'VariableNames', {'P_Pa', 'P_psi', 'phases', 'beta_vapour', 'V_over_Vsat', ...
                                  'SL', 'rhoM_L', 'rhoM_V', 'reason'});
        end

        function [tbl, compositions] = cvd(obj, T, z, r_cap, Pdew, Plist)
            % Classical CVD on a dew-point fluid: at each stage flash, then
            % remove vapour at constant P until the cell is back at V_cell.
            % Pdew [Pa]: saturation (dew) pressure; Plist [Pa] < Pdew.
            % compositions: struct array of remaining z, liquid x, produced y.
            arguments
                obj
                T     (1,1) double {mustBePositive}
                z     (:,1) double {mustBeNonnegative}
                r_cap (1,1) double {mustBePositive}
                Pdew  (1,1) double {mustBePositive}
                Plist (:,1) double {mustBePositive}
            end
            z = z / sum(z);
            [~, ~, V0] = obj.EOS.calculateState(Pdew, T, z, -1, r_cap);
            Vcell = V0;                     % 1 mol of feed at Pdew
            nTot  = 1.0;
            zCur  = z;

            Plist = sort(Plist(Plist < Pdew), 'descend');
            n = numel(Plist);
            [SL, gasProd, nLeft, beta] = deal(nan(n, 1));
            phases = zeros(n, 1);
            reason = strings(n, 1);
            compositions = repmat(struct('P', NaN, 'z', z, 'x', z, 'y', z), n, 1);
            cumGas = 0.0;
            K = [];
            for i = 1:n
                r = obj.Flash.solve(Plist(i), T, zCur, r_cap, 'K_seed', K);
                phases(i) = r.phases;  reason(i) = r.reason;  beta(i) = r.beta;
                if ~r.converged
                    warning('PVTExperiments:FlashFailed', ...
                        'CVD stopped at %.1f psi: flash %s.', Plist(i) / 6894.757, r.reason);
                    break;                                      % later stages depend on this one
                end
                if r.phases == 1
                    % Single phase at this stage: expand, remove excess as feed-composition fluid
                    Vt = nTot * r.V_total;
                    dn = max(Vt - Vcell, 0) / r.V_total;
                    nTot = nTot - dn;
                    cumGas = cumGas + dn;
                    SL(i) = double(r.label == "liquid");
                    compositions(i) = struct('P', Plist(i), 'z', zCur, 'x', zCur, 'y', zCur);
                else
                    K  = r.K;
                    nL = nTot * (1 - r.beta);
                    nV = nTot * r.beta;
                    VL = nL * r.V_L;
                    VV = nV * r.V_V;
                    SL(i) = VL / Vcell;
                    dn = max(VL + VV - Vcell, 0) / r.V_V;      % vapour moles removed
                    if dn > nV
                        error('PVTExperiments:CellOverfilled', ...
                              'Liquid volume exceeds the cell volume at %.1f psi.', Plist(i) / 6894.757);
                    end
                    nV = nV - dn;
                    nTot = nL + nV;
                    cumGas = cumGas + dn;
                    compositions(i) = struct('P', Plist(i), 'z', zCur, 'x', r.x, 'y', r.y);
                    zCur = (nL * r.x + nV * r.y) / nTot;
                end
                gasProd(i) = cumGas;
                nLeft(i) = nTot;
            end
            tbl = table(Plist, Plist / 6894.757, phases, beta, SL, gasProd, nLeft, reason, ...
                'VariableNames', {'P_Pa', 'P_psi', 'phases', 'beta_vapour', 'SL_cell', ...
                                  'cumGasProduced_molFrac', 'molesLeft', 'reason'});
        end
    end
end
