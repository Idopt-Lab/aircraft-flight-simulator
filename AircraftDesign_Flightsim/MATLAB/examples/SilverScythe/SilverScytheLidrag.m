function out = SilverScytheLidrag(eta,ccl_ca,AR,n_harmonics)
% SILVERSCYTHELIDRAG  Lift-induced drag and span efficiency from a spanload.
%
% MATLAB implementation of the Trefftz-plane analysis performed by Dave
% Ives' and Dr. William H. Mason's LIDRAG program, which the report used to
% obtain the Silver Scythe span efficiency (report 11.4).  LIDRAG expands a
% symmetric spanload in odd Fourier (Glauert/Munk) harmonics and evaluates
% the induced drag in the Trefftz plane:
%
%   eta = cos(theta),  ccl/ca = 4*AR * sum_{n odd} A_n sin(n theta)
%   CL   = pi*AR*A_1
%   CD_i = pi*AR*sum_{n odd} n*A_n^2
%   e    = CL^2/(pi*AR*CD_i) = A_1^2 / sum_{n odd} n*A_n^2
%
% A spanload with non-zero tip loading has a formally divergent harmonic
% series, so LIDRAG - like every discrete implementation - regularises the
% problem through a finite harmonic count.  n_harmonics odd terms
% (n = 1,3,...,2*n_harmonics-1) are fitted to the supplied stations in a
% least-squares sense.  The default of 9 terms was chosen by reproducing
% four independent LIDRAG runs of the executable:
%
%   spanload            LIDRAG      this function (9 terms)
%   elliptic            1.00000     1.00000
%   Silver Scythe       0.89127     0.88925   (-0.23%)
%   triangular          0.71988     0.72190   (+0.28%)
%   rectangular         0.60504     0.60072   (-0.71%)
%
% The archived Silver Scythe LIDRAG output is kept in
% reference/lidrag_silverscythe.out and is checked by SilverScytheSelfCheck.
%
% Inputs
%   eta          spanwise stations y/(b/2), root to tip, in [0,1]
%   ccl_ca       spanload c*cl/c_a at those stations
%   AR           wing aspect ratio
%   n_harmonics  number of odd harmonics (default 9)
%
% Output struct
%   CL, CD_i, e_span, k_induced, A (harmonic amplitudes), n (harmonic orders)

if nargin < 4 || isempty(n_harmonics), n_harmonics = 9; end

eta = eta(:);
ccl_ca = ccl_ca(:);
if numel(eta) ~= numel(ccl_ca)
    error('SilverScytheLidrag:SizeMismatch', ...
        'eta and ccl_ca must have the same length.');
end
if numel(eta) < 2*n_harmonics
    n_harmonics = max(floor(numel(eta)/2),1);
end
if ~(isscalar(AR) && isfinite(AR) && AR > 0)
    error('SilverScytheLidrag:InvalidAR','AR must be a positive scalar.');
end

theta = acos(min(max(eta,0),1));
n = (1:2:(2*n_harmonics-1));

% ccl/ca = 4*AR*sum A_n sin(n theta); solve for B_n = 4*AR*A_n
Basis = sin(theta*n);
B = Basis\ccl_ca;

denominator = sum(n(:).*B(:).^2);
if denominator <= 0
    out = struct('CL',0,'CD_i',0,'e_span',1,'k_induced',1/(pi*AR), ...
        'A',zeros(numel(n),1),'n',n(:),'method',"LIDRAG (Trefftz plane)");
    return;
end

CL     = (pi/4)*B(1);
e_span = B(1)^2/denominator;
e_span = min(max(e_span,1e-3),1.0);
k      = 1/(pi*AR*e_span);

out = struct();
out.CL        = CL;
out.e_span    = e_span;
out.k_induced = k;
out.CD_i      = k*CL^2;
out.A         = B(:)/(4*AR);
out.n         = n(:);
out.n_harmonics = n_harmonics;
out.method    = "LIDRAG (Ives/Mason) Trefftz-plane harmonic analysis";

end
