function out = SilverScytheThrustModel(thr,mach,alt_m,V_mps,mode,ss,eng)
% SILVERSCYTHETHRUSTMODEL  Single-engine F414-GE-400 thrust and fuel flow.
%
% Signature is the one TurbofanPropulsion calls: (throttle, Mach,
% altitude_m, airspeed).  mode, ss and eng are bound by
% SilverScytheContext, so the handle handed to the propulsion element is
%
%   @(thr,M,alt,V) SilverScytheThrustModel(thr,M,alt,V,mode,ss,eng)
%
% Throttle convention, matching the report's engine decks:
%   0.0  idle
%   0.5  maximum dry (military)
%   1.0  maximum afterburning
%
% mode = "mattingly"  (low fidelity)
%   The report's own engine deck, section 13.5: Mattingly's algebraic
%   thrust lapse in theta0 and delta0, with the F414-GE-400 sea-level
%   static values of 13,000 lbf dry and 22,000 lbf wet and TSFC of 0.84 and
%   1.85.  This is a per-engine version of the senior-design
%   computeThrust414.m, which used the twin-engine totals.
%
% mode = "npss"  (high fidelity)
%   The team's archived NPSS deck (report section 13.6) wherever it has
%   coverage, which is Mach <= 0.60 and altitude <= 35,000 ft.  Outside
%   that box the Mattingly lapse is used, multiplied by the NPSS/Mattingly
%   thrust and TSFC ratios evaluated at the nearest point on the coverage
%   boundary.  The result is continuous across the boundary and equals NPSS
%   exactly inside it, so takeoff, catapult launch, approach and low-level
%   climb - every condition the report's carrier-suitability requirements
%   are written against - run on measured-cycle data, while the supersonic
%   dash falls back to the validated algebraic lapse.
%
% Output struct consumed by TurbofanPropulsion
%   thrust      [N]        fuel_flow  [kg/s]
%   direction   local unit thrust axis
%   moment      local intrinsic moment
% plus diagnostics: thrust_lbf, TSFC, pla, source.

if nargin < 4 || isempty(V_mps), V_mps = 0; end %#ok<NASGU>
if nargin < 5 || isempty(mode),  mode = "mattingly"; end
if nargin < 6 || isempty(ss),    ss = SilverScytheData(); end
mode = lower(string(mode));

thr = max(0,min(1,thr));
mach = max(0,mach);
alt_ft = max(alt_m,0)/ss.units.ft_to_m;

switch mode
    case "mattingly"
        [T_lbf,TSFC] = mattingly_point(ss,thr,mach,alt_ft);
        pla = throttle_to_pla(thr,10,99,150);
        source = "Mattingly algebraic deck (report 13.5)";

    case "npss"
        if nargin < 7 || isempty(eng)
            eng = SilverScytheEngineDeck(ss);
        end
        pla = throttle_to_pla(thr,eng.pla_idle,eng.pla_mil,eng.pla_max_ab);

        m_in = min(mach,eng.mach_max);
        a_in = min(max(alt_ft,eng.alt_ft_grid(1)),eng.alt_max_ft);

        T_npss  = eng.Fn_fn(m_in,a_in,pla);
        TSFC_npss = eng.TSFC_fn(m_in,a_in,pla);

        inside = (mach <= eng.mach_max + 1e-9) && ...
                 (alt_ft <= eng.alt_max_ft + 1e-9);

        if inside
            T_lbf = T_npss;
            TSFC  = TSFC_npss;
            source = "NPSS deck (report 13.6)";
        else
            [T_edge,TSFC_edge] = mattingly_point(ss,thr,m_in,a_in);
            [T_here,TSFC_here] = mattingly_point(ss,thr,mach,alt_ft);
            kT = T_npss/max(T_edge,1e-6);
            kS = TSFC_npss/max(TSFC_edge,1e-6);
            % Keep the boundary correction bounded so a single odd cell
            % cannot distort the supersonic end of the envelope.
            kT = min(max(kT,0.5),2.0);
            kS = min(max(kS,0.5),2.0);
            T_lbf = T_here*kT;
            TSFC  = TSFC_here*kS;
            source = "Mattingly scaled to NPSS boundary";
        end

    otherwise
        error('SilverScytheThrustModel:UnknownMode', ...
            'mode must be "mattingly" or "npss".');
end

T_lbf = max(T_lbf,0);
TSFC  = max(TSFC,0);

out = struct();
out.thrust     = T_lbf*ss.units.lbf_to_N;
out.fuel_flow  = TSFC*T_lbf*ss.units.lbm_to_kg/3600;   % kg/s
out.direction  = [1;0;0];
out.moment     = zeros(3,1);
out.thrust_lbf = T_lbf;
out.TSFC       = TSFC;
out.pla        = pla;
out.mode       = mode;
out.source     = source;

end

% =======================================================================
function pla = throttle_to_pla(thr,pla_idle,pla_mil,pla_ab)
if thr <= 0.5
    pla = pla_idle + (thr/0.5)*(pla_mil - pla_idle);
else
    pla = pla_mil + ((thr - 0.5)/0.5)*(pla_ab - pla_mil);
end
end

% =======================================================================
function [T_lbf,TSFC] = mattingly_point(ss,thr,mach,alt_ft)
% Per-engine version of computeThrust414.m.

gamma = ss.engine.mattingly.gamma;
TR    = ss.engine.mattingly.TR;

[T_sl_K,~,P_sl_Pa,~] = atmosisa(0);
Tstd = T_sl_K*1.8;
Pstd = P_sl_Pa*0.0208854342;

[T_K,~,P_Pa,~] = atmosisa(alt_ft*ss.units.ft_to_m);
Ti = T_K*1.8;
Pi = P_Pa*0.0208854342;

theta  = Ti/Tstd;
delta  = Pi/Pstd;
theta0 = theta*(1 + ((gamma-1)*mach^2)/2);
delta0 = delta*(1 + ((gamma-1)*mach^2)/2)^(gamma/(gamma-1));

if theta0 <= TR
    lapse_dry = delta0*(1 - 0.3*mach);
    lapse_wet = delta0*(1 - 0.1*sqrt(mach));
else
    lapse_dry = delta0*(1 - 0.3*mach - 1.7*(theta0-TR)/theta0);
    lapse_wet = delta0*(1 - 0.1*sqrt(mach) - 2.2*(theta0-TR)/theta0);
end

T_dry = max(ss.engine.T_SL_dry_lbf*lapse_dry,0);
T_wet = max(ss.engine.T_SL_AB_lbf *lapse_wet,0);

TSFC_dry = ss.engine.TSFC_SL_dry*(1 + 0.35*abs(mach))*sqrt(theta);
TSFC_wet = ss.engine.TSFC_SL_AB *(1 + 0.35*abs(mach - 0.4))*sqrt(theta);

% Idle taken as 8% of the dry rating; the report's deck has no idle line.
T_idle = 0.08*T_dry;

if thr <= 0.5
    s = thr/0.5;
    T_lbf = T_idle + s*(T_dry - T_idle);
    TSFC  = TSFC_dry;
else
    s = (thr - 0.5)/0.5;
    T_lbf = T_dry + s*(T_wet - T_dry);
    TSFC  = TSFC_dry + s*(TSFC_wet - TSFC_dry);
end

end
