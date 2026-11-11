function out = SilverScytheSpanload(ss,n_stations)
% SILVERSCYTHESPANLOAD  Normalised spanload fed to LIDRAG.
%
% MATLAB port of the spanload construction in Aero/SUMAC/lift_calc.m and of
% report section 11.4: the wing is divided into 50 spanwise stations, the
% section lift coefficient is taken as constant across the span (the report
% multiplies the XFOIL 2-D lift slope by the effective angle of attack),
% and the local load is normalised by the average chord,
%
%   ccl/ca = c(eta)*cl / c_a,     c_a = S_ref/b
%
% Because the section cl is constant, the *shape* of the spanload - and
% therefore the LIDRAG span efficiency - depends only on the chord
% distribution.  That is why the Silver Scythe span efficiency is a single
% planform constant rather than a function of flight condition, exactly as
% the archived LIDRAG output shows.
%
% Output struct
%   eta, chord_ft, ccl_ca (for unit section cl), c_a_ft

if nargin < 1 || isempty(ss), ss = SilverScytheData(); end
if nargin < 2 || isempty(n_stations), n_stations = ss.deck.spanload_stations; end

eta = linspace(0,1,n_stations).';
chord_ft = (ss.tip_chord_ft - ss.root_chord_ft)*eta + ss.root_chord_ft;
c_a_ft = ss.mean_chord_ft;

out = struct();
out.eta      = eta;
out.chord_ft = chord_ft;
out.c_a_ft   = c_a_ft;
out.ccl_ca   = chord_ft/c_a_ft;   % unit section cl
out.method   = "lift_calc.m spanload, 50 stations, constant section cl";

end
