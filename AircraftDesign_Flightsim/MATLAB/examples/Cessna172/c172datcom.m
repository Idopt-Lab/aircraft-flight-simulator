function lookup = c172datcom(datcom_file_path)
if nargin < 1 || isempty(datcom_file_path)
    error('DATCOM file path required');
end
if ~exist(datcom_file_path, 'file')
    error('File not found: %s', datcom_file_path);
end

data = parse_C172datcom(datcom_file_path);
lookup = @(x, u, geom) eval_C172datcom(x, u, geom, data);
end

function coeff = eval_C172datcom(x, u, geom, d) %#ok<INUSD>
u_b = double(x(4));
v_b = double(x(5));
w_b = double(x(6));

V = sqrt(u_b^2 + v_b^2 + w_b^2);
alpha = atan2(w_b, max(abs(u_b), 1e-9));
beta = asin(max(-1, min(1, v_b / max(V, 1e-9))));

u_d = double(u(:));
da = 0;
de = 0;
dr = 0;

if numel(u_d) >= 1, da = u_d(1); end
if numel(u_d) >= 2, de = u_d(2); end
if numel(u_d) >= 3, dr = u_d(3); end

CL_base = interp1(d.alpha, d.CL, alpha, 'linear', 'extrap');
CD_base = interp1(d.alpha, d.CD, alpha, 'linear', 'extrap');
CM_base = interp1(d.alpha, d.CM, alpha, 'linear', 'extrap');

coeff = struct();
coeff.CL = CL_base + d.CLDE * de;
coeff.CD = CD_base + d.CD0_add + d.CD_alpha_extra * alpha^2;
coeff.CY = d.CYB * beta + d.CYDR * dr;
coeff.Cl = d.CLB * beta + d.CLDA * da;
coeff.Cm = CM_base + d.CMDE * de;
coeff.Cn = d.CNB * beta + d.CNDA * da + d.CNDR * dr;

coeff.CL = max(min(coeff.CL, d.CL_max_abs), -d.CL_max_abs);

coeff.CLA = d.CLA;
coeff.CMA = d.CMA;
coeff.CYB = d.CYB;
coeff.CNB = d.CNB;
coeff.CLB = d.CLB;
coeff.CLQ = d.CLQ;
coeff.CMQ = d.CMQ;
coeff.CLP = d.CLP;
coeff.CYP = d.CYP;
coeff.CNP = d.CNP;
coeff.CNR = d.CNR;
coeff.CLR = d.CLR;

coeff.CLDE = d.CLDE;
coeff.CMDE = d.CMDE;
coeff.CLDA = d.CLDA;
coeff.CNDA = d.CNDA;
coeff.CYDR = d.CYDR;
coeff.CNDR = d.CNDR;

coeff.takeoff_params = d.takeoff_params;
coeff.landing_params = d.landing_params;
end

function data = parse_C172datcom(filepath)
text = fileread(filepath);
lines = strsplit(text, '\n');

alpha_vec = [];
CL_vec = [];
CD_vec = [];
CM_vec = [];

CLA = NaN;
CMA = NaN;
CYB = NaN;
CNB = NaN;
CLB = NaN;
CLQ = NaN;
CMQ = NaN;
CLP = NaN;
CNR = NaN;
CLR = NaN;
CYP = NaN;
CNP = NaN;

full_config_tag = 'WING-BODY-VERTICAL TAIL-HORIZONTAL TAIL CONFIGURATION';
damp_tag = 'DYNAMIC DERIVATIVES (PER RADIAN)';
in_full_config = false;
full_config_done = false;
in_damp_full = false;
damp_full_done = false;

i = 1;
while i <= numel(lines)
    line = lines{i};

    if contains(line, full_config_tag)
        in_full_config = true;
        in_damp_full = false;
    end

    if in_full_config && ~full_config_done
        if contains(line, 'ALPHA') && contains(line, 'CD') && contains(line, 'CL') && contains(line, 'CM')
            k = i + 1;
            while k <= numel(lines) && isempty(strtrim(lines{k}))
                k = k + 1;
            end

            temp_alpha = [];
            temp_CL = [];
            temp_CD = [];
            temp_CM = [];
            temp_CLA = NaN;
            temp_CMA = NaN;
            temp_CYB = NaN;
            temp_CNB = NaN;
            temp_CLB = NaN;

            while k <= numel(lines)
                dataline = strtrim(lines{k});

                if isempty(dataline)
                    k = k + 1;
                    continue;
                end

                if strncmp(dataline, '1', 1) || strncmp(dataline, '0***', 4)
                    break;
                end

                vals = str2num(dataline); %#ok<ST2NM>
                if ~isempty(vals) && numel(vals) >= 4
                    a = vals(1);
                    CD = vals(2);
                    CL = vals(3);
                    CM = vals(4);

                    if abs(a) < 30
                        temp_alpha(end+1,1) = a; %#ok<AGROW>
                        temp_CL(end+1,1) = CL; %#ok<AGROW>
                        temp_CD(end+1,1) = CD; %#ok<AGROW>
                        temp_CM(end+1,1) = CM; %#ok<AGROW>

                        if numel(vals) >= 8  && abs(a) <= 2.5 && abs(vals(8))  > 0.01  && isnan(temp_CLA), temp_CLA = vals(8);  end
                        if numel(vals) >= 9  && abs(a) <= 2.5 && abs(vals(9))  > 0.01  && isnan(temp_CMA), temp_CMA = vals(9);  end
                        if numel(vals) >= 10 && abs(a) <= 2.5 && abs(vals(10)) > 0.001 && isnan(temp_CYB), temp_CYB = vals(10); end
                        if numel(vals) >= 11 && abs(a) <= 2.5 && abs(vals(11)) > 0.001 && isnan(temp_CNB), temp_CNB = vals(11); end
                        if numel(vals) >= 12 && isnan(temp_CLB), temp_CLB = vals(12); end
                    end
                end
                k = k + 1;
            end

            if numel(temp_alpha) > numel(alpha_vec)
                alpha_vec = temp_alpha;
                CL_vec = temp_CL;
                CD_vec = temp_CD;
                CM_vec = temp_CM;

                if ~isnan(temp_CLA), CLA = temp_CLA; end
                if ~isnan(temp_CMA), CMA = temp_CMA; end
                if ~isnan(temp_CYB), CYB = temp_CYB; end
                if ~isnan(temp_CNB), CNB = temp_CNB; end
                if ~isnan(temp_CLB), CLB = temp_CLB; end

                full_config_done = true;
            end
            in_full_config = false;
        end
    end

    if contains(line, full_config_tag) && full_config_done
        in_damp_full = true;
    end

    if in_damp_full && ~damp_full_done && contains(line, damp_tag)
        for j = i+4:min(i+20, numel(lines))
            parts = strsplit(strtrim(lines{j}));
            parts = parts(~cellfun('isempty', parts));

            if numel(parts) >= 10
                a = str2double(parts{1});
                if isfinite(a) && abs(a) <= 1
                    if isnan(CLQ), CLQ = str2double(parts{2}); end
                    if isnan(CMQ), CMQ = str2double(parts{3}); end
                    if isnan(CLP), CLP = str2double(parts{6}); end
                    if isnan(CYP), CYP = str2double(parts{7}); end
                    if isnan(CNP), CNP = str2double(parts{8}); end
                    if isnan(CNR), CNR = str2double(parts{9}); end
                    if isnan(CLR), CLR = str2double(parts{10}); end

                    damp_full_done = true;
                    break;
                end
            end
        end
        in_damp_full = false;
    end

    i = i + 1;
end

if isnan(CLA), CLA = 5.2; end
if isnan(CMA), CMA = -0.8; end
if isnan(CYB), CYB = -0.35; end
if isnan(CNB), CNB = 0.08; end
if isnan(CLB), CLB = -0.10; end
if isnan(CLQ), CLQ = 3.5; end
if isnan(CMQ), CMQ = -10.0; end
if isnan(CLP), CLP = -0.50; end
if isnan(CNR), CNR = -0.10; end
if isnan(CLR), CLR = 0.15; end
if isnan(CYP), CYP = -0.20; end
if isnan(CNP), CNP = -0.04; end

[a_u, idx] = unique(alpha_vec, 'sorted');
data = struct();
data.alpha = deg2rad(a_u(:));
data.CL = CL_vec(idx);
data.CD = CD_vec(idx);
data.CM = CM_vec(idx);

data.CLA = CLA;
data.CMA = CMA;
data.CYB = CYB;
data.CNB = CNB;
data.CLB = CLB;
data.CLQ = CLQ;
data.CMQ = CMQ;
data.CLP = CLP;
data.CYP = CYP;
data.CNP = CNP;
data.CNR = CNR;
data.CLR = CLR;

data.CLDE = 0.30;
data.CMDE = -1.20;
data.CLDA = 0.075;
data.CNDA = 0.008;
data.CYDR = 0.20;
data.CNDR = -0.075;

data.CD0_add = 0.0;
data.CD_alpha_extra = 0.0;
data.CL_max_abs = 1.9;

data.takeoff_params = struct( ...
    'mu_rolling', 0.03, ...
    'mu_braking', 0.35, ...
    'safety_factor', 1.15, ...
    'CLmax_takeoff', 1.6, ...
    'VR_to_Vs_ratio', 1.10, ...
    'V2_to_Vs_ratio', 1.20, ...
    'screen_height_takeoff_ft', 50, ...
    'rotation_alpha_deg', 8, ...
    'reaction_time_s', 1.0);

data.landing_params = struct( ...
    'mu_braking', 0.50, ...
    'approach_angle_deg', 3.0, ...
    'safety_factor', 1.67, ...
    'CLmax_landing', 1.9, ...
    'Vapp_to_Vs_ratio', 1.30, ...
    'screen_height_landing_ft', 50, ...
    'flare_height_m', 3.0, ...
    'idle_throttle', 0.05, ...
    'use_idle_thrust', true);
end