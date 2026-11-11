function D = SilverScytheDrag(ctx,CL,mach,altitude_m)
% SILVERSCYTHEDRAG  Silver Scythe drag polar for both fidelity levels.
%
% low   The textbook drag polar of the senior-design script
%       drag_calc_sol.m (Raymer / Roskam / Nikolai):
%         CD_min   = Cfe*S_wet/S_ref                       (Raymer 12.3)
%         e_osw    = 4.6(1-0.033 AR^0.53) cos(LE)^0.1 - 3.3
%         k1_sub   = 1/(pi e_osw AR)
%         CL_alpha = a0/(1 + 57.3 a0/(pi e_alpha AR))
%         CL_minD  = CL_alpha*(-alpha_L0)/2
%         k2_sub   = -2 k1_sub CL_minD
%         k1_sup   = AR(M^2-1)/(4 AR sqrt(M^2-1) - 2) cos(LE)
%         CD_wave  = 4.5 pi/S (Amax/l)^2 EWD (0.74+0.37 cos LE)
%                    (1 - 0.3 sqrt(M - M_CD0))
%       with the same 0.95-1.05 transonic blend.
%
% high  The report's higher-fidelity chain, tabulated by
%       SilverScytheAeroDeck:
%         CD_profile  Mason's FRICTION, per-component, Raymer form and
%                     interference factors, 5% stealth leakage
%         CD_induced  k*CL^2 with k = 1/(pi AR e), e from Ives'/Mason's
%                     LIDRAG applied to the 50-station spanload
%         CD_wave     Raymer drag-divergence and Sears-Haack method
%
% Inputs
%   ctx         SilverScytheContext struct
%   CL          lift coefficient, scalar or vector
%   mach        Mach number, scalar or matching CL
%   altitude_m  altitude [m], scalar or matching CL
%
% Output struct
%   CD_profile, CD_induced, CD_wave, CD_total, k_induced, e_span, CL_minD

ss = ctx.ss;

CL = CL(:).';
if isscalar(mach), mach = mach*ones(size(CL)); end
if isscalar(altitude_m), altitude_m = altitude_m*ones(size(CL)); end
mach = mach(:).';
altitude_m = altitude_m(:).';

alt_ft = altitude_m/ss.units.ft_to_m;
n = numel(CL);

CD_profile = zeros(1,n);
CD_induced = zeros(1,n);
CD_wave    = zeros(1,n);
k_out      = zeros(1,n);
e_out      = zeros(1,n);

switch ctx.fidelity
    case "low"
        p = ctx.drag_low;
        CL_minD_out = zeros(1,n);
        for i = 1:n
            [k1,CLminD,CDmin,CDw] = low_coefficients(p,mach(i));
            CD_profile(i) = CDmin;
            CD_induced(i) = k1*(CL(i) - CLminD)^2;
            CD_wave(i)    = CDw;
            k_out(i)      = k1;
            e_out(i)      = 1/(pi*ss.AR*max(k1,1e-12));
            CL_minD_out(i) = CLminD;
        end
        CL_minD = CL_minD_out;
        source = "textbook drag polar (drag_calc_sol.m: Raymer/Roskam/Nikolai)";

    case "high"
        deck = ctx.deck;
        m_q = min(max(mach,deck.mach_grid(1)),deck.mach_grid(end));
        a_q = min(max(alt_ft,deck.alt_ft_grid(1)),deck.alt_ft_grid(end));
        CL_q = min(max(CL,deck.CL_grid(1)),deck.CL_grid(end));

        CD_profile = deck.CD_profile_fn(m_q(:),a_q(:)).';
        CD_wave    = deck.CD_wave_fn(m_q(:),CL_q(:)).';
        CD_induced = max(deck.k_induced*CL.^2,0);
        k_out(:)   = deck.k_induced;
        e_out(:)   = deck.e_span;
        CL_minD    = zeros(1,n);
        source = "FRICTION + LIDRAG + Raymer wave drag (report 11.4-11.5)";

    otherwise
        error('SilverScytheDrag:UnknownFidelity', ...
            'ctx.fidelity must be "low" or "high".');
end

CD_total = max(CD_profile + CD_induced + CD_wave,1e-4);

D = struct();
D.CL         = CL;
D.mach       = mach;
D.altitude_m = altitude_m;
D.CD_profile = CD_profile;
D.CD_induced = CD_induced;
D.CD_wave    = CD_wave;
D.CD_total   = CD_total;
D.k_induced  = k_out;
D.e_span     = e_out;
D.CL_minD    = CL_minD;
D.source     = source;

end

% -----------------------------------------------------------------------
function [k1,CLminD,CDmin,CDwave] = low_coefficients(p,M)
% Textbook subsonic/supersonic coefficients with the drag_calc_sol.m blend.
%
% drag_calc_sol.m writes the polar as CD = CD0 + k1*CL^2 + k2*CL with
% CD0 = CDmin + k1*CL_minD^2 and k2 = -2*k1*CL_minD.  That is the expanded
% form of CD = CDmin + k1*(CL - CL_minD)^2, which is used here because the
% completed square cannot go negative for any CL.  The CL_minD offset is
% dropped supersonically, matching k2_sup = 0 in the original script.

k1_sub = p.k1_sub;
CDmin_sub = p.CDmin_sub;

Msup = max(M,1.0001);
denominator = 4*p.AR*sqrt(Msup^2 - 1) - 2;
k1_sup = max((p.AR*(Msup^2 - 1)/max(denominator,1e-9))*cosd(p.LE_sweep_deg),0);
CLminD_sup = 0;

if M < p.M_cutoff
    wave_term = p.wave_multiplier*(1 - 0.3*sqrt(max(M - p.M_CD0,0)));
else
    wave_term = p.wave_multiplier*(1 - 0.3*sqrt(max(p.M_cutoff - p.M_CD0,0)));
end
if M <= p.M_CD0
    wave_term = 0;
end
wave_term = max(wave_term,0);
CDmin_sup = CDmin_sub;

M1 = p.blend(1);
M2 = p.blend(2);

if M < M1
    k1 = k1_sub;  CLminD = p.CL_minD;  CDmin = CDmin_sub;  CDwave = 0;
elseif M > M2
    k1 = k1_sup;  CLminD = CLminD_sup; CDmin = CDmin_sup;  CDwave = wave_term;
else
    w = (M - M1)/(M2 - M1);
    k1     = (1-w)*k1_sub    + w*k1_sup;
    CLminD = (1-w)*p.CL_minD + w*CLminD_sup;
    CDmin  = (1-w)*CDmin_sub + w*CDmin_sup;
    CDwave = w*wave_term;
end

end
