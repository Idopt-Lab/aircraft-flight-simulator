function lookup = example_aero_parser(filepath)
    if ~exist(filepath,'file'), error('File not found'); end
    data = parse_file(filepath);
    lookup = @(x, u, geom) evaluate(x, u, geom, data);
end

function data = parse_file(filepath)
    text = fileread(filepath);
    lines = strsplit(text, '\n');
    
    alpha = [];
    CL = [];
    CD = [];
    CM = [];
    
    d = struct('CLA',NaN,'CMA',NaN,'CYB',NaN,'CNB',NaN,'CLB',NaN,...
               'CLQ',NaN,'CMQ',NaN,'CLP',NaN,'CNR',NaN,'CLR',NaN,'CYP',NaN,'CNP',NaN);
    
    for i = 1:length(lines)
        line = lines{i};
        
        if contains(line, 'DATA_TABLE_MARKER')
            for j = i:min(i+100, length(lines))
                if contains(lines{j}, 'ALPHA') && contains(lines{j}, 'CL')
                    k = j + 1;
                    
                    while k <= length(lines)
                        vals = str2num(strtrim(lines{k}));
                        
                        if isempty(vals) || contains(lines{k}, 'END')
                            break;
                        end
                        
                        if length(vals) >= 4 && abs(vals(1)) < 30
                            alpha(end+1) = vals(1);
                            CD(end+1) = vals(2);
                            CL(end+1) = vals(3);
                            CM(end+1) = vals(4);
                            
                            if abs(vals(1)) <= 2.5 && length(vals) >= 9
                                if isnan(d.CLA) && vals(8) > 2, d.CLA = vals(8); end
                                if isnan(d.CMA), d.CMA = vals(9); end
                            end
                        end
                        k = k + 1;
                    end
                    break;
                end
            end
        end
        
        if contains(line, 'DYNAMIC')
            for j = i+4:min(i+10, length(lines))
                parts = strsplit(strtrim(lines{j}));
                parts(cellfun('isempty',parts)) = [];
                
                if length(parts) >= 10
                    a = str2double(parts{1});
                    if isfinite(a) && abs(a) <= 1
                        if isnan(d.CLQ), d.CLQ = str2double(parts{2}); end
                        if isnan(d.CMQ), d.CMQ = str2double(parts{3}); end
                        if isnan(d.CLP), d.CLP = str2double(parts{6}); end
                        if isnan(d.CNR), d.CNR = str2double(parts{9}); end
                        break;
                    end
                end
            end
        end
    end
    
    defaults = struct('CLA',5.0,'CMA',-0.8,'CYB',-0.3,'CNB',0.02,'CLB',-0.19,...
                     'CLQ',7.3,'CMQ',-14.7,'CLP',-0.46,'CNR',-0.05,'CLR',0.04,'CYP',-0.05,'CNP',-0.01);
    
    fn = fieldnames(defaults);
    for k = 1:length(fn)
        if isnan(d.(fn{k})), d.(fn{k}) = defaults.(fn{k}); end
    end
    
    [au, ix] = unique(alpha);
    data.alpha = deg2rad(au(:));
    data.CL = CL(ix)';
    data.CD = CD(ix)';
    data.CM = CM(ix)';
    
    fn = fieldnames(d);
    for k = 1:length(fn)
        data.(fn{k}) = d.(fn{k});
    end
    
    data.CLDE = 0.30;
    data.CMDE = -1.20;
    data.CLDA = 0.075;
    data.CNDA = 0.008;
    data.CYDR = 0.20;
    data.CNDR = -0.075;
end

function C = evaluate(x, u, geom, data)
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
    
    CL = interp1(data.alpha, data.CL, alpha, 'linear', 'extrap') + data.CLDE * de;
    CD = interp1(data.alpha, data.CD, alpha, 'linear', 'extrap');
    Cm = interp1(data.alpha, data.CM, alpha, 'linear', 'extrap') + data.CMDE * de;
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