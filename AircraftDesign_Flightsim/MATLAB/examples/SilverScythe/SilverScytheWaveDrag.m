function out = SilverScytheWaveDrag(ss,mach,CL,EWD)
% SILVERSCYTHEWAVEDRAG  Transonic and supersonic wave drag.
%
% MATLAB port of the senior-design Aero/SUMAC/emp_wave.m, which implements
% the report's method (section 11.5): a Raymer curve fit for the zero-lift
% drag-divergence Mach number (report Figure 21), a Raymer lift correction
% (report Figure 22), historical curve fits through the transonic drag rise
% between M_cr and Mach 1.2, and a Sears-Haack equivalent-body estimate
% above Mach 1.2 (report Equations 15 and 16).
%
% Two implementation notes, both preserved deliberately so that this
% function reproduces the published numbers:
%
%  1. emp_wave.m evaluates the Sears-Haack term as (9*pi/2)*(Amax/l), i.e.
%     with the area ratio to the first power, whereas report Equation 15 as
%     typeset squares it.  The un-squared form is what produced the
%     published Figure 24 value (0.016106 against 0.0161335 at Mach 1.64),
%     so it is the form used here.
%  2. The transonic curve fits are high-order polynomials with severe
%     cancellation and can return small negative values near Mach 0.95.
%     The result is clamped at zero, which is the only change made to the
%     original expressions.
%
% Inputs
%   ss     SilverScytheData struct
%   mach   Mach number, scalar or vector
%   CL     lift coefficient, scalar or the same size as mach
%   EWD    empirical wave-drag factor (default ss.wave.EWD = 1.8)
%
% Output struct
%   CD_wave, M_DD, M_cr, quarter_chord_sweep_deg, regime

if nargin < 1 || isempty(ss), ss = SilverScytheData(); end
if nargin < 4 || isempty(EWD), EWD = ss.wave.EWD; end

mach = mach(:).';
if isscalar(CL), CL = CL*ones(size(mach)); end
CL = CL(:).';
if numel(CL) ~= numel(mach)
    error('SilverScytheWaveDrag:SizeMismatch','CL must match mach in size.');
end

if abs(ss.tc_wing - 0.06) > 1e-9
    error('SilverScytheWaveDrag:UnsupportedThickness', ...
        'The Raymer curve fits ported here are tabulated for t/c = 0.06 only.');
end

b     = ss.span_ft;
cr    = ss.root_chord_ft;
ct    = ss.tip_chord_ft;
S_ref = ss.S_ref_ft2;
LamLE = ss.LE_sweep_deg;
Amax  = ss.max_cross_section_ft2;
l_ref = ss.wave_reference_length_ft;

% Quarter-chord sweep, emp_wave.m expression
Lambda_c4 = atand(tand(LamLE) - (cr - ct)/(2*b));

% Zero-lift drag-divergence Mach number, Raymer fit for t/c = 0.06
M_DD_0L = -0.000000043366233*Lambda_c4^4 ...
          +0.000006652171359*Lambda_c4^3 ...
          -0.000336136773529*Lambda_c4^2 ...
          +0.009329824226811*Lambda_c4 ...
          +0.737692981629277;

n = numel(mach);
CD_wave = zeros(1,n);
M_DD    = zeros(1,n);
M_cr    = zeros(1,n);
regime  = strings(1,n);

for i = 1:n
    M = mach(i);
    CLi = CL(i);

    % Raymer lift correction factor, t/c = 0.06 fit
    LFDD =  0.195466937024053*CLi^4 ...
          + 0.037170857674685*CLi^3 ...
          - 0.246230576285308*CLi^2 ...
          + 0.014721829243993*CLi ...
          + 0.998673837852488;

    M_DD(i) = M_DD_0L*LFDD - 0.05*CLi;
    M_cr(i) = M_DD(i) - 0.08;

    if M > 1.2
        Dq_SH   = (9*pi)/2*(Amax/l_ref);
        Dq_wave = EWD*(1 - 0.386*(M - 1.2)^0.57* ...
                  (1 - (pi*LamLE^0.77)/100))*Dq_SH;
        CD_wave(i) = Dq_wave/S_ref;
        regime(i)  = "supersonic Sears-Haack";
    elseif M > 1.05
        CD_wave(i) =  2.787903576788589*M^3 ...
                    - 9.644552898645681*M^2 ...
                    +11.103597445886850*M ...
                    - 4.227956762127869;
        regime(i) = "curve fit 1.05-1.20";
    elseif M > 1.0
        CD_wave(i) = 0.242010299205642*M - 0.229399215061641;
        regime(i)  = "linear fit 1.00-1.05";
    elseif M > M_cr(i)
        CD_wave(i) =  25.626857418690157*M^4 ...
                    - 89.174087773670720*M^3 ...
                    +116.467261762068300*M^2 ...
                    - 67.646485325736847*M ...
                    + 14.739024005808796;
        regime(i) = "curve fit M_cr-1.00";
    else
        CD_wave(i) = 0;
        regime(i)  = "subcritical";
    end
end

CD_wave(imag(CD_wave) ~= 0) = 0;
CD_wave = max(real(CD_wave),0);

out = struct();
out.mach    = mach;
out.CL      = CL;
out.CD_wave = CD_wave;
out.M_DD    = M_DD;
out.M_cr    = M_cr;
out.M_DD_zero_lift = M_DD_0L;
out.quarter_chord_sweep_deg = Lambda_c4;
out.EWD     = EWD;
out.regime  = regime;
out.method  = "emp_wave.m (Raymer/Brandt), report Eqs 15-16";

end
