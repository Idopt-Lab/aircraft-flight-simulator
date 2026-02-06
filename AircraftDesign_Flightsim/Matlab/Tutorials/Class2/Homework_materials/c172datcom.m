function lookup = c172datcom(datcom_file_path)
    if nargin < 1 || isempty(datcom_file_path)
        error('c172datcom:NoFile','DATCOM file path required');
    end
    if ~exist(datcom_file_path,'file')
        error('c172datcom:FileNotFound','File not found: %s', datcom_file_path);
    end

    data = parse_datcom(datcom_file_path);
    lookup = @(state_vec, control_vec, geometry) eval_datcom(state_vec, control_vec, geometry, data);
end

function data = parse_datcom(filepath)
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

    i = 1;
    while i <= length(lines)
        line = lines{i};

        if contains(line, 'ALPHA') && contains(line, 'CD') && contains(line, 'CL') && contains(line, 'CM')
            k = i + 1;
            while k <= length(lines) && isempty(strtrim(lines{k}))
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

            while k <= length(lines)
                dataline = strtrim(lines{k});

                if isempty(dataline)
                    k = k + 1;
                    continue;
                end

                if startsWith(dataline, '1') || startsWith(dataline, '0***')
                    break;
                end

                vals = str2num(dataline);

                if ~isempty(vals) && length(vals) >= 4
                    alpha = vals(1);
                    CD = vals(2);
                    CL = vals(3);
                    CM = vals(4);

                    if abs(alpha) < 30
                        temp_alpha(end+1) = alpha;
                        temp_CL(end+1) = CL;
                        temp_CD(end+1) = CD;
                        temp_CM(end+1) = CM;

                        if length(vals) >= 8 && abs(alpha) <= 2.5 && abs(vals(8)) > 0.01 && isnan(temp_CLA)
                            temp_CLA = vals(8);
                        end

                        if length(vals) >= 9 && abs(alpha) <= 2.5 && abs(vals(9)) > 0.01 && isnan(temp_CMA)
                            temp_CMA = vals(9);
                        end

                        if length(vals) >= 10 && abs(alpha) <= 2.5 && abs(vals(10)) > 0.001 && isnan(temp_CYB)
                            temp_CYB = vals(10);
                        end

                        if length(vals) >= 11 && abs(alpha) <= 2.5 && abs(vals(11)) > 0.001 && isnan(temp_CNB)
                            temp_CNB = vals(11);
                        end

                        if length(vals) >= 12 && isnan(temp_CLB)
                            temp_CLB = vals(12);
                        end
                    end
                end

                k = k + 1;
            end

            if ~isempty(temp_alpha) && length(temp_alpha) > length(alpha_vec)
                alpha_vec = temp_alpha;
                CL_vec = temp_CL;
                CD_vec = temp_CD;
                CM_vec = temp_CM;
                if ~isnan(temp_CLA), CLA = temp_CLA; end
                if ~isnan(temp_CMA), CMA = temp_CMA; end
                if ~isnan(temp_CYB), CYB = temp_CYB; end
                if ~isnan(temp_CNB), CNB = temp_CNB; end
                if ~isnan(temp_CLB), CLB = temp_CLB; end
            end
        end

        if contains(line, 'DYNAMIC DERIVATIVES (PER RADIAN)')
            for j = i+4:min(i+20, length(lines))
                dynline = strtrim(lines{j});
                parts = strsplit(dynline);
                parts = parts(~cellfun('isempty', parts));

                if length(parts) >= 10
                    alpha_dyn = str2double(parts{1});
                    if isfinite(alpha_dyn) && abs(alpha_dyn) <= 1.0
                        if isnan(CLQ), CLQ = str2double(parts{2}); end
                        if isnan(CMQ), CMQ = str2double(parts{3}); end
                        if isnan(CLP), CLP = str2double(parts{6}); end
                        if isnan(CYP), CYP = str2double(parts{7}); end
                        if isnan(CNP), CNP = str2double(parts{8}); end
                        if isnan(CNR), CNR = str2double(parts{9}); end
                        if isnan(CLR), CLR = str2double(parts{10}); end
                        break;
                    end
                end
            end
        end

        i = i + 1;
    end

    if isempty(alpha_vec)
        error('parse_datcom:NoData','No aerodynamic data found');
    end

    if isnan(CLA), CLA = 4.7; end
    if isnan(CMA), CMA = -0.9; end
    if isnan(CYB), CYB = -0.31; end
    if isnan(CNB), CNB = 0.071; end
    if isnan(CLB), CLB = -0.089; end
    if isnan(CLQ), CLQ = 3.9; end
    if isnan(CMQ), CMQ = -12.4; end
    if isnan(CLP), CLP = -0.49; end
    if isnan(CNR), CNR = -0.095; end
    if isnan(CLR), CLR = 0.14; end
    if isnan(CYP), CYP = -0.19; end
    if isnan(CNP), CNP = -0.03; end

    [alpha_unique, idx] = unique(alpha_vec);

    data = struct();
    data.alpha = deg2rad(alpha_unique(:));
    data.CL = CL_vec(idx)';
    data.CD = CD_vec(idx)';
    data.CM = CM_vec(idx)';

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
end

function C = eval_datcom(x, u, geom, data)
    vel = x(4:6);
    omega = x(10:12);

    V = max(norm(vel), 1e-6);
    alpha = atan2(vel(3), max(abs(vel(1)), 1e-9));
    beta = atan2(vel(2), max(sqrt(vel(1)^2 + vel(3)^2), 1e-9));

    de = 0; if length(u) >= 2, de = u(2); end
    da = 0; if length(u) >= 1, da = u(1); end
    dr = 0; if length(u) >= 3, dr = u(3); end

    b = geom.wing_span;
    c = geom.mean_aerodynamic_chord;

    CL_base = interp1(data.alpha, data.CL, alpha, 'linear', 'extrap');
    CD_base = interp1(data.alpha, data.CD, alpha, 'linear', 'extrap');
    Cm_base = interp1(data.alpha, data.CM, alpha, 'linear', 'extrap');

    CL = CL_base + data.CLDE * de;
    CD = CD_base;
    Cm = Cm_base + data.CMDE * de;

    CY = data.CYB * beta + data.CYDR * dr;
    Cl = data.CLB * beta + data.CLDA * da;
    Cn = data.CNB * beta + data.CNDR * dr + data.CNDA * da;

    if V > 1
        p_hat = omega(1) * b / (2 * V);
        q_hat = omega(2) * c / (2 * V);
        r_hat = omega(3) * b / (2 * V);

        CL = CL + data.CLQ * q_hat;
        Cm = Cm + data.CMQ * q_hat;
        Cl = Cl + data.CLP * p_hat + data.CLR * r_hat;
        Cn = Cn + data.CNP * p_hat + data.CNR * r_hat;
        CY = CY + data.CYP * p_hat;
    end

    C = struct('CL', CL, 'CD', max(CD, 0), 'CY', CY, 'Cl', Cl, 'Cm', Cm, 'Cn', Cn);
end