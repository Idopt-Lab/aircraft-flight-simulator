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
        
        if contains(line, 'CASEID') && contains(line, 'Vertical Tail') && contains(line, 'Horizontal Tail')
            
            for j = i:min(i+200, length(lines))
                header = lines{j};
                
                if contains(header, 'ALPHA') && contains(header, 'CL') && contains(header, 'CD')
                    
                    k = j + 1;
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
                        
                        if startsWith(dataline, '1') || contains(dataline, 'FLIGHT') || contains(dataline, 'DYNAMIC')
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
                                
                                if length(vals) >= 8 && abs(alpha) <= 2.5 && vals(8) > 2.0 && isnan(temp_CLA)
                                    temp_CLA = vals(8);
                                end
                                
                                if length(vals) >= 9 && abs(alpha) <= 2.5 && abs(vals(9)) > 0.05 && isnan(temp_CMA)
                                    temp_CMA = vals(9);
                                end
                                
                                if length(vals) >= 10 && abs(alpha) <= 2.5 && isnan(temp_CYB)
                                    temp_CYB = vals(10);
                                end
                                
                                if length(vals) >= 11 && abs(alpha) <= 2.5 && isnan(temp_CNB)
                                    temp_CNB = vals(11);
                                end
                                
                                if length(vals) >= 12 && abs(alpha) <= 2.5 && isnan(temp_CLB)
                                    temp_CLB = vals(12);
                                end
                            end
                        end
                        
                        k = k + 1;
                    end
                    
                    if ~isempty(temp_alpha)
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
                    
                    break;
                end
            end
        end
        
        if contains(line, 'DYNAMIC DERIVATIVES (PER RADIAN)')
            for j = i+4:min(i+10, length(lines))
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
    
    if isnan(CLA), CLA = 5.0; fprintf('Warning: CLA missing, using default 5.0\n'); end
    if isnan(CMA), CMA = -0.8; fprintf('Warning: CMA missing, using default -0.8\n'); end
    if isnan(CYB), CYB = -0.3; fprintf('Warning: CYB missing, using default -0.3\n'); end
    if isnan(CNB), CNB = 0.02; fprintf('Warning: CNB missing, using default 0.02\n'); end
    if isnan(CLB), CLB = -0.19; fprintf('Warning: CLB missing, using default -0.19\n'); end
    if isnan(CLQ), CLQ = 7.3; fprintf('Warning: CLQ missing, using default 7.3\n'); end
    if isnan(CMQ), CMQ = -14.7; fprintf('Warning: CMQ missing, using default -14.7\n'); end
    if isnan(CLP), CLP = -0.46; fprintf('Warning: CLP missing, using default -0.46\n'); end
    if isnan(CNR), CNR = -0.05; fprintf('Warning: CNR missing, using default -0.05\n'); end
    if isnan(CLR), CLR = 0.04; fprintf('Warning: CLR missing, using default 0.04\n'); end
    if isnan(CYP), CYP = -0.05; fprintf('Warning: CYP missing, using default -0.05\n'); end
    if isnan(CNP), CNP = -0.01; fprintf('Warning: CNP missing, using default -0.01\n'); end
    
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
    
    fprintf('\n=== DATCOM Data Loaded ===\n');
    fprintf('Alpha: %.1f to %.1f deg (%d points)\n', rad2deg(min(data.alpha)), rad2deg(max(data.alpha)), length(data.alpha));
    fprintf('CLA=%.3f  CMA=%.3f\n', data.CLA, data.CMA);
    fprintf('CYB=%.3f  CNB=%.3f  CLB=%.3f\n', data.CYB, data.CNB, data.CLB);
    fprintf('CLQ=%.3f  CMQ=%.3f\n', data.CLQ, data.CMQ);
    fprintf('CLP=%.3f  CNR=%.3f  CLR=%.3f\n', data.CLP, data.CNR, data.CLR);
    fprintf('CYP=%.3f  CNP=%.3f\n', data.CYP, data.CNP);
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