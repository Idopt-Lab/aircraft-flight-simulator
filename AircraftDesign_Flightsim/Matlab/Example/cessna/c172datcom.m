function lookup = c172datcom(datcom_file_path)
    if nargin < 1 || isempty(datcom_file_path)
        error('DATCOM file path required');
    end
    if ~exist(datcom_file_path,'file')
        error('File not found');
    end
    data = parse_datcom(datcom_file_path);
    lookup = @(state_vec, control_vec, geometry) eval_datcom(state_vec, control_vec, geometry, data);
end

function data = parse_datcom(filepath)
    text = fileread(filepath);
    lines = strsplit(text, '\n');

    TARGET_CONFIG = 'WING-BODY-HORIZONTAL TAIL';

    alpha_vec = []; CL_vec = []; CD_vec = []; CM_vec = [];
    CLA = NaN; CMA = NaN; CYB = NaN; CNB = NaN; CLB = NaN;
    CLQ = NaN; CMQ = NaN; CLP = NaN; CNR = NaN; CLR = NaN; CYP = NaN; CNP = NaN;

    i = 1;
    while i <= numel(lines)
        line = lines{i};

        if contains(line,'ALPHA') && contains(line,'CD') && contains(line,'CL') && contains(line,'CM')
            config_label = '';
            for back = max(1,i-15):i-1
                if contains(lines{back},'CONFIGURATION') || contains(lines{back},'WING') || contains(lines{back},'BODY')
                    config_label = strtrim(lines{back});
                    break;
                end
            end

            is_target = contains(config_label, TARGET_CONFIG) && ...
                        ~contains(config_label, 'VERTICAL') && ...
                        ~startsWith(config_label, 'BODY') && ...
                        ~startsWith(config_label, 'WING ALONE') && ...
                        ~startsWith(config_label, 'HORIZONTAL TAIL');

            if ~is_target
                i = i + 1;
                continue;
            end

            k = i + 1;
            while k <= numel(lines) && isempty(strtrim(lines{k}))
                k = k + 1;
            end

            temp_alpha = []; temp_CL = []; temp_CD = []; temp_CM = [];
            temp_CLA = NaN; temp_CMA = NaN; temp_CYB = NaN; temp_CNB = NaN; temp_CLB = NaN;

            while k <= numel(lines)
                dataline = strtrim(lines{k});
                if isempty(dataline), k = k+1; continue; end
                if startsWith(dataline,'1') || startsWith(dataline,'0***'), break; end

                vals = str2num(dataline);
                if ~isempty(vals) && numel(vals) >= 4
                    alp = vals(1); CD = vals(2); CL = vals(3); CM = vals(4);
                    if abs(alp) < 30
                        temp_alpha(end+1) = alp;
                        temp_CL(end+1)    = CL;
                        temp_CD(end+1)    = CD;
                        temp_CM(end+1)    = CM;
                        if numel(vals)>=8  && abs(alp)<=4.0 && abs(vals(8)) >0.01  && isnan(temp_CLA), temp_CLA=vals(8);  end
                        if numel(vals)>=9  && abs(alp)<=4.0 && abs(vals(9)) >0.01  && isnan(temp_CMA), temp_CMA=vals(9);  end
                        if numel(vals)>=10 && abs(alp)<=4.0 && abs(vals(10))>0.001 && isnan(temp_CYB), temp_CYB=vals(10); end
                        if numel(vals)>=11 && abs(alp)<=4.0 && abs(vals(11))>0.001 && isnan(temp_CNB), temp_CNB=vals(11); end
                        if numel(vals)>=12 && isnan(temp_CLB), temp_CLB=vals(12); end
                    end
                end
                k = k + 1;
            end

            if ~isempty(temp_alpha)
                alpha_vec = temp_alpha; CL_vec = temp_CL; CD_vec = temp_CD; CM_vec = temp_CM;
                if ~isnan(temp_CLA), CLA=temp_CLA; end
                if ~isnan(temp_CMA), CMA=temp_CMA; end
                if ~isnan(temp_CYB), CYB=temp_CYB; end
                if ~isnan(temp_CNB), CNB=temp_CNB; end
                if ~isnan(temp_CLB), CLB=temp_CLB; end
            end
        end

        if contains(line,'DYNAMIC DERIVATIVES (PER RADIAN)')
            for j = i+4:min(i+20,numel(lines))
                parts = strsplit(strtrim(lines{j}));
                parts = parts(~cellfun('isempty',parts));
                if numel(parts) >= 10
                    a = str2double(parts{1});
                    if isfinite(a) && abs(a) <= 1
                        if isnan(CLQ), CLQ=str2double(parts{2}); end
                        if isnan(CMQ), CMQ=str2double(parts{3}); end
                        if isnan(CLP), CLP=str2double(parts{6}); end
                        if isnan(CYP), CYP=str2double(parts{7}); end
                        if isnan(CNP), CNP=str2double(parts{8}); end
                        if isnan(CNR), CNR=str2double(parts{9}); end
                        if isnan(CLR), CLR=str2double(parts{10}); end
                        break;
                    end
                end
            end
        end
        i = i + 1;
    end

    if isempty(alpha_vec)
        error('c172datcom: WING-BODY-HORIZONTAL TAIL table not found in %s', filepath);
    end

    if isnan(CLA), CLA= 5.2;   end
    if isnan(CMA), CMA=-0.8;   end
    if isnan(CYB), CYB=-0.35;  end
    if isnan(CNB), CNB= 0.08;  end
    if isnan(CLB), CLB=-0.10;  end
    if isnan(CLQ), CLQ= 3.5;   end
    if isnan(CMQ), CMQ=-10.0;  end
    if isnan(CLP), CLP=-0.50;  end
    if isnan(CNR), CNR=-0.10;  end
    if isnan(CLR), CLR= 0.15;  end
    if isnan(CYP), CYP=-0.20;  end
    if isnan(CNP), CNP=-0.04;  end

    [a_u, idx] = unique(alpha_vec);
    data.alpha = deg2rad(a_u(:));
    data.CL    = CL_vec(idx)';
    data.CD    = CD_vec(idx)';
    data.CM    = CM_vec(idx)';

    data.CLA=CLA; data.CMA=CMA; data.CYB=CYB; data.CNB=CNB; data.CLB=CLB;
    data.CLQ=CLQ; data.CMQ=CMQ; data.CLP=CLP; data.CYP=CYP; data.CNP=CNP; data.CNR=CNR; data.CLR=CLR;

    data.CLDE           = 0.30;
    data.CMDE           = -1.20;
    data.CLDA           = 0.075;
    data.CNDA           = 0.008;
    data.CYDR           = 0.20;
    data.CNDR           = -0.075;
    data.CD0_add        = 0.055;
    data.CD_alpha_extra = 0.004;
    data.CL_max_abs     = 1.9;
    data.alpha_max_rad  = deg2rad(25);
end

function C = eval_datcom(x, u, geom, data)
    vel   = x(4:6);
    omega = x(10:12);
    V     = max(norm(vel), 1e-6);
    u_b   = vel(1); v_b = vel(2); w_b = vel(3);

    alpha = atan2(w_b, max(u_b, 1e-9));
    beta  = atan2(v_b, max(sqrt(u_b^2+w_b^2), 1e-9));

    de = 0; if numel(u)>=2, de=u(2); end
    da = 0; if numel(u)>=1, da=u(1); end
    dr = 0; if numel(u)>=3, dr=u(3); end

    c = geom.mean_aerodynamic_chord;
    b = geom.wing_span;

    a  = min(max(alpha, -data.alpha_max_rad), data.alpha_max_rad);
    CL = interp1(data.alpha, data.CL, a, 'linear', 'extrap') + data.CLDE*de;
    CD = interp1(data.alpha, data.CD, a, 'linear', 'extrap') + data.CD0_add + data.CD_alpha_extra*abs(a);
    Cm = interp1(data.alpha, data.CM, a, 'linear', 'extrap') + data.CMDE*de;

    CL = max(min(CL, data.CL_max_abs), -data.CL_max_abs);
    CD = max(CD, 0.01);

    CY = data.CYB*beta  + data.CYDR*dr;
    Cl = data.CLB*beta  + data.CLDA*da;
    Cn = data.CNB*beta  + data.CNDR*dr + data.CNDA*da;

    if V > 1
        p_hat = omega(1)*b/(2*V);
        q_hat = omega(2)*c/(2*V);
        r_hat = omega(3)*b/(2*V);
        CL = CL + data.CLQ*q_hat;
        Cm = Cm + data.CMQ*q_hat;
        Cl = Cl + data.CLP*p_hat + data.CLR*r_hat;
        Cn = Cn + data.CNP*p_hat + data.CNR*r_hat;
        CY = CY + data.CYP*p_hat;
    end

    takeoff_params = struct('mu_rolling',0.03,'mu_braking',0.35,'safety_factor',1.15, ...
        'CLmax_takeoff',1.6,'VR_to_Vs_ratio',1.10,'V2_to_Vs_ratio',1.20, ...
        'screen_height_takeoff_ft',50,'rotation_alpha_deg',8,'reaction_time_s',1.0);
    landing_params = struct('mu_braking',0.50,'approach_angle_deg',3.0,'safety_factor',1.67, ...
        'CLmax_landing',1.9,'Vapp_to_Vs_ratio',1.30,'screen_height_landing_ft',50, ...
        'flare_height_m',3.0,'idle_throttle',0.05,'use_idle_thrust',true);

    C = struct('CL',CL,'CD',CD,'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn, ...
               'takeoff_params',takeoff_params,'landing_params',landing_params);
end