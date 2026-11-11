function airfoil = SilverScytheAirfoil(ss)
% SILVERSCYTHEAIRFOIL  NACA 64-206 section data from the XFOIL polar.
%
% MATLAB port of the senior-design Aero/cl_alpha.m.  The 2-D lift-curve
% slope and zero-lift angle are obtained from a first-order fit through the
% complete XFOIL polar, exactly as cl_alpha.m does, so the values feeding
% the high-fidelity lift model are the ones the report's analyses used.
%
% Output fields
%   alpha_deg, CL, CD, CM   raw XFOIL polar, sorted, duplicates removed
%   cl_alpha_per_deg        2-D lift slope from the full-range fit  [1/deg]
%   cl_alpha_per_rad        same, per radian
%   alpha_L0_deg            zero-lift angle of attack              [deg]
%   cl_at_zero_alpha        section CL at alpha = 0
%   alpha_max_deg           clean stall angle used by the team
%   cl_max_2D               maximum section CL present in the polar
%   linear_fit              slope/alpha_L0 from the -4..8 deg window only
%
% Reference: Aero/cl_alpha.m, Gobbler Ghost Program 2026.

if nargin < 1 || isempty(ss)
    ss = SilverScytheData();
end

persistent cache cache_key
key = char(fullfile(ss.paths.data_dir,char(ss.airfoil.polar_file)));
if ~isempty(cache) && strcmp(cache_key,key)
    airfoil = cache;
    return;
end

if ~isfile(key)
    error('SilverScytheAirfoil:MissingPolar', ...
        'XFOIL polar not found: %s',key);
end

A = readmatrix(key,'FileType','text');
A = A(all(isfinite(A(:,1:5)),2),1:5);      % alpha CL CD CDp CM
A = sortrows(A,1);
[~,idx] = unique(A(:,1));
A = A(idx,:);

alpha = A(:,1);
cl    = A(:,2);

% Full-range first-order fit: this is what cl_alpha.m uses.
fit_full = polyfit(alpha,cl,1);

% Linear-window fit, reported for comparison only.
window = alpha >= -4 & alpha <= 8;
fit_lin = polyfit(alpha(window),cl(window),1);

airfoil = struct();
airfoil.name             = ss.airfoil.name;
airfoil.polar_file       = key;
airfoil.alpha_deg        = alpha;
airfoil.CL               = cl;
airfoil.CD               = A(:,3);
airfoil.CM               = A(:,5);
airfoil.cl_alpha_per_deg = fit_full(1);
airfoil.cl_alpha_per_rad = fit_full(1)*180/pi;
airfoil.alpha_L0_deg     = -fit_full(2)/fit_full(1);
airfoil.cl_at_zero_alpha = fit_full(2);
airfoil.alpha_min_deg    = min(alpha);
airfoil.alpha_max_deg    = ss.airfoil.alpha_max_deg;
airfoil.cl_max_2D        = max(cl);
airfoil.linear_fit       = struct( ...
    'cl_alpha_per_deg',fit_lin(1), ...
    'alpha_L0_deg',-fit_lin(2)/fit_lin(1));

cache = airfoil;
cache_key = key;

end
