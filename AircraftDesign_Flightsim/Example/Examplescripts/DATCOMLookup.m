function lookup = DATCOMLookup(datcom_file_path)
if nargin < 1 || isempty(datcom_file_path)
    error('DATCOMLookup:FileRequired', 'DATCOM file path required.');
end
if ~exist(datcom_file_path,'file')
    error('DATCOMLookup:NotFound', 'DATCOM file not found: %s', datcom_file_path);
end
data = read_datcom(datcom_file_path);
lookup = @(state_vec, control_vec, geometry) datcom_eval(state_vec, control_vec, geometry, data);
end
function data = read_datcom(fp)
if exist('datcomimport','file') ~= 2
    error('DATCOMLookup:MissingToolbox', ...
        'datcomimport not found. Install/enable MATLAB Aerospace Toolbox (or add your own parser).');
end
raw = datcomimport(fp);
if ~iscell(raw), raw = {raw}; end
if isempty(raw), error('DATCOMLookup:Empty', 'No cases found in DATCOM file.'); end
c = raw{1};
deg2rad_ = pi/180;
data = struct();
data.alpha = [];
data.CLtab = [];
data.CDtab = [];
data.CMtab = [];
if isfield(c,'alpha') && isnumeric(c.alpha)
    data.alpha = c.alpha(:) * deg2rad_;
end
if isfield(c,'cl') && isnumeric(c.cl), data.CLtab = c.cl(:); end
if isfield(c,'cd') && isnumeric(c.cd), data.CDtab = c.cd(:); end
if isfield(c,'cm') && isnumeric(c.cm), data.CMtab = c.cm(:); end
data.der = struct();
if isfield(c,'cla'), data.der.CLA = c.cla; end
if isfield(c,'cma'), data.der.CMA = c.cma; end
if isfield(c,'cyb'), data.der.CYB = c.cyb; end
if isfield(c,'clb'), data.der.CLB = c.clb; end
if isfield(c,'cnb'), data.der.CNB = c.cnb; end
data.ctrl = struct();
if isfield(c,'clroll'), data.ctrl.CLDA = c.clroll / deg2rad_; end
if isfield(c,'cn_asy'), data.ctrl.CNDA = c.cn_asy / deg2rad_; end
data.CL0 = interp0(data.alpha, data.CLtab, 0.10);
data.CD0 = interp0(data.alpha, data.CDtab, 0.025);
data.CM0 = interp0(data.alpha, data.CMtab, 0.02);
if ~isfield(data.der,'CLA'), data.der.CLA = 4.5; end
if ~isfield(data.der,'CMA'), data.der.CMA = -0.6; end
if ~isfield(data.der,'CYB'), data.der.CYB = -0.5; end
if ~isfield(data.der,'CLB'), data.der.CLB = -0.12; end
if ~isfield(data.der,'CNB'), data.der.CNB = 0.12; end
end
function C = datcom_eval(x, u, geom, data)
vel   = x(4:6);
omega = x(10:12);
V = max(norm(vel), 1e-6);
alpha = atan2(vel(3), max(vel(1),1e-9));
beta  = atan2(vel(2), max(sqrt(vel(1)^2 + vel(3)^2),1e-9));
da = get_u(u,1);
de = get_u(u,2);
dr = get_u(u,3);
S = geom.wing_area;
b = geom.wing_span;
if isprop(geom,'mean_aerodynamic_chord')
    c = geom.mean_aerodynamic_chord;
elseif isfield(geom,'ref_chord')
    c = geom.ref_chord;
else
    c = max(S/max(b,1e-6), 1e-6);
end
CL = data.CL0 + data.der.CLA*alpha;
CD = data.CD0;
Cm = data.CM0 + data.der.CMA*alpha;
CY = data.der.CYB*beta;
Cl = data.der.CLB*beta;
Cn = data.der.CNB*beta;
if isfield(data.ctrl,'CLDA'), Cl = Cl + data.ctrl.CLDA*da; end
if isfield(data.ctrl,'CNDA'), Cn = Cn + data.ctrl.CNDA*da; end
% If you have elevator/rudder derivatives from DATCOM, plug them here:
% CL = CL + CLDE*de;  Cm = Cm + CMDE*de;  CY = CY + CYDR*dr; etc.
AR = b^2 / max(S,1e-6);
e  = 0.80;
K  = 1/(pi*e*AR);
CD = CD + K*CL^2;
p_hat = omega(1)*b/(2*V);
q_hat = omega(2)*c/(2*V);
r_hat = omega(3)*b/(2*V);
CLq = 3.0;
Cmq = -10.0;
Clp = -0.50;
Cnr = -0.20;
CL = CL + CLq*q_hat;
Cm = Cm + Cmq*q_hat;
Cl = Cl + Clp*p_hat;
Cn = Cn + Cnr*r_hat;
takeoff_params = struct('mu_ground',0.03,'CLmax_takeoff',1.6,'CD0_takeoff',0.040);
landing_params = struct('mu_braking',0.35,'CLmax_landing',2.0,'CD0_landing',0.055);
C = struct('CL',CL,'CD',max(CD,0),'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn, ...
           'takeoff_params',takeoff_params,'landing_params',landing_params);

end
function y0 = interp0(x, y, default_val)
if isempty(x) || isempty(y) || numel(x) ~= numel(y)
    y0 = default_val;
    return;
end
y0 = interp1(x, y, 0, 'linear', 'extrap');
end

function v = get_u(u, idx)
v = 0;
if isnumeric(u) && numel(u) >= idx
    v = u(idx);
end
end
