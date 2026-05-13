function lookup = b737datcom(datcom_file_path)
if nargin < 1 || isempty(datcom_file_path)
    error('DATCOM file path required');
end
if ~exist(datcom_file_path, 'file')
    error('File not found: %s', datcom_file_path);
end

txt   = fileread(datcom_file_path);
txt   = strrep(txt, sprintf('\r\n'), sprintf('\n'));
txt   = strrep(txt, sprintf('\r'),   sprintf('\n'));
lines = strsplit(txt, '\n');

ft_to_m = 0.3048;

data.control = struct( ...
    'CLDE',0.00, 'CMDE',-1.10, ...
    'CYDR',0.18, 'CNDR',-0.070, ...
    'CLDA',0.070,'CNDA',0.010);

data.alpha_max_rad = deg2rad(25);
data.CL_max_abs    = 2.2;
data.CD_min        = 0.01;

data.takeoff_params = struct( ...
    'mu_rolling',0.03,'mu_braking',0.35,'safety_factor',1.15, ...
    'CLmax_takeoff',1.8,'VR_to_Vs_ratio',1.10,'V2_to_Vs_ratio',1.20, ...
    'screen_height_takeoff_ft',35,'rotation_alpha_deg',8,'reaction_time_s',1.0);

data.landing_params = struct( ...
    'mu_braking',0.45,'approach_angle_deg',3.0,'safety_factor',1.67, ...
    'CLmax_landing',2.2,'Vapp_to_Vs_ratio',1.30,'screen_height_landing_ft',50, ...
    'flare_height_m',3.0,'idle_throttle',0.05,'use_idle_thrust',true);

case_defs = extract_case_definitions(lines, ft_to_m);
tables    = extract_aero_tables(lines);

if isempty(tables)
    error('No aerodynamic tables found in DATCOM output.');
end

fprintf('\n=== DATCOM TABLES EXTRACTED ===\n');
for i = 1:numel(tables)
    fprintf('%2d) M=%.3f  Alt=%.0f ft  points=%d  config=%s\n', ...
        i, tables(i).mach, tables(i).altitude_ft, ...
        numel(tables(i).alpha_deg), char(tables(i).config_label));
end

data.grid = build_aero_grid(tables, case_defs, ft_to_m);

fprintf('\nAerodynamic grid: %d alpha x %d Mach points\n', ...
    numel(data.grid.alpha_rad), numel(data.grid.mach_vec));
fprintf('Mach range:  %.3f - %.3f\n', min(data.grid.mach_vec), max(data.grid.mach_vec));
fprintf('Alpha range: %.1f - %.1f deg\n', ...
    rad2deg(min(data.grid.alpha_rad)), rad2deg(max(data.grid.alpha_rad)));

lookup = @(x, u, geom) local_eval_lookup(x, u, geom, data);
end

% ── Build unified alpha-Mach grid ─────────────────────────────────────────────

function grid = build_aero_grid(tables, case_defs, ft_to_m) %#ok<INUSD>

if numel(tables) == 1
    m0    = tables(1).mach;
    beta0 = sqrt(max(1 - m0^2, 0.01));
    mach_extra = [0.20 0.30 0.40 0.50 0.60 0.70 0.74 0.78];
    mach_extra = mach_extra(abs(mach_extra - m0) > 0.01);
    orig = tables(1);
    for im = 1:numel(mach_extra)
        me     = mach_extra(im);
        beta_e = sqrt(max(1 - me^2, 0.01));
        t      = orig;
        t.mach = me;
        t.CL   = orig.CL * (beta0 / beta_e);
        t.CM   = orig.CM * (beta0 / beta_e);
        Mdd    = 0.74;
        cd_wave  = 0;
        if me > Mdd, cd_wave = 20*(me-Mdd)^4; end
        cd_press = max(orig.CD - 0.008, 0);
        cd_skin  = min(orig.CD, 0.008);
        t.CD   = cd_skin + cd_press*(beta0/beta_e) + cd_wave;
        t.config_priority = orig.config_priority - 0.1;
        tables(end+1) = t; %#ok<AGROW>
    end
    fprintf('  (single table — Prandtl-Glauert scaled to %d Mach points)\n', numel(tables));
end

mach_vec         = unique(round([tables.mach], 4));
alpha_deg_common = unique(vertcat(tables.alpha_deg));
alpha_rad        = deg2rad(alpha_deg_common);
n_a = numel(alpha_rad);
n_m = numel(mach_vec);

CL_mat = NaN(n_a, n_m);
CD_mat = NaN(n_a, n_m);
CM_mat = NaN(n_a, n_m);

for im = 1:n_m
    m          = mach_vec(im);
    candidates = tables(abs([tables.mach] - m) < 0.005);
    if isempty(candidates), continue; end
    t = candidates(end);
    for ia = 1:n_a
        a = alpha_deg_common(ia);
        if a >= min(t.alpha_deg) && a <= max(t.alpha_deg)
            CL_mat(ia,im) = interp1(t.alpha_deg, t.CL, a, 'linear');
            CD_mat(ia,im) = interp1(t.alpha_deg, t.CD, a, 'linear');
            CM_mat(ia,im) = interp1(t.alpha_deg, t.CM, a, 'linear');
        end
    end
end

for im = 1:n_m
    valid = ~isnan(CL_mat(:,im));
    if sum(valid) >= 2
        av = alpha_rad(valid);
        CL_mat(:,im) = interp1(av, CL_mat(valid,im), alpha_rad, 'linear', 'extrap');
        CD_mat(:,im) = interp1(av, CD_mat(valid,im), alpha_rad, 'linear', 'extrap');
        CM_mat(:,im) = interp1(av, CM_mat(valid,im), alpha_rad, 'linear', 'extrap');
    end
end

grid.alpha_rad = alpha_rad;
grid.mach_vec  = mach_vec;
grid.CL        = CL_mat;
grid.CD        = CD_mat;
grid.CM        = CM_mat;

grid.CYB = -0.45; grid.CNB =  0.09; grid.CLB = -0.08;
grid.CLQ =  3.5;  grid.CMQ = -12.0; grid.CLP = -0.55;
grid.CYP = -0.20; grid.CNP = -0.05; grid.CNR = -0.12; grid.CLR =  0.10;
end

% ── Runtime lookup ─────────────────────────────────────────────────────────────

function C = local_eval_lookup(x, u, geom, data)
vel   = x(4:6);
omega = x(10:12);
V     = max(norm(vel), 1e-6);
u_b   = vel(1); v_b = vel(2); w_b = vel(3);

alpha = atan2(w_b, max(u_b, 1e-9));
beta  = atan2(v_b, max(sqrt(u_b^2 + w_b^2), 1e-9));

alt_m = max(-x(3), 0);
[~, a_spd, ~, ~] = atmosisa(alt_m);
mach  = V / max(a_spd, 1e-9);

g    = data.grid;
b    = max(geom.wing_span,              1e-9);
cbar = max(geom.mean_aerodynamic_chord, 1e-9);

a_clip = min(max(alpha, -data.alpha_max_rad), data.alpha_max_rad);
m_clip = min(max(mach, min(g.mach_vec)), max(g.mach_vec));

if numel(g.mach_vec) >= 2
    CL = interp2(g.mach_vec, g.alpha_rad, g.CL, m_clip, a_clip, 'linear');
    CD = interp2(g.mach_vec, g.alpha_rad, g.CD, m_clip, a_clip, 'linear');
    Cm = interp2(g.mach_vec, g.alpha_rad, g.CM, m_clip, a_clip, 'linear');
else
    CL = interp1(g.alpha_rad, g.CL(:,1), a_clip, 'linear', 'extrap');
    CD = interp1(g.alpha_rad, g.CD(:,1), a_clip, 'linear', 'extrap');
    Cm = interp1(g.alpha_rad, g.CM(:,1), a_clip, 'linear', 'extrap');
end

if ~isfinite(CL)
    CL = interp1(g.alpha_rad, g.CL(:,end), a_clip, 'linear', 'extrap');
    CD = interp1(g.alpha_rad, g.CD(:,end), a_clip, 'linear', 'extrap');
    Cm = interp1(g.alpha_rad, g.CM(:,end), a_clip, 'linear', 'extrap');
end

da = 0; de = 0; dr = 0;
if numel(u) >= 1, da = u(1); end
if numel(u) >= 2, de = u(2); end
if numel(u) >= 3, dr = u(3); end

CL = CL + data.control.CLDE * de;
Cm = Cm + data.control.CMDE * de;
CY = g.CYB * beta + data.control.CYDR * dr;
Cl = g.CLB * beta + data.control.CLDA * da;
Cn = g.CNB * beta + data.control.CNDR * dr + data.control.CNDA * da;

if V > 1
    p_hat = omega(1)*b    / (2*V);
    q_hat = omega(2)*cbar / (2*V);
    r_hat = omega(3)*b    / (2*V);
    CL = CL + g.CLQ*q_hat;
    Cm = Cm + g.CMQ*q_hat;
    Cl = Cl + g.CLP*p_hat + g.CLR*r_hat;
    Cn = Cn + g.CNP*p_hat + g.CNR*r_hat;
    CY = CY + g.CYP*p_hat;
end

CL = max(min(CL, data.CL_max_abs), -data.CL_max_abs);
CD = max(CD, data.CD_min);

C = struct('CL',CL,'CD',CD,'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn, ...
    'takeoff_params',data.takeoff_params,'landing_params',data.landing_params);
end

% ── Case definition extraction ─────────────────────────────────────────────────

function case_defs = extract_case_definitions(lines, ft_to_m)
case_defs = struct('caseid',{},'mach',{},'altitude_ft',{},'gamma_deg',{}, ...
    'thrust_lbf_per_engine',{},'phase',{},'altitude_m',{},'gamma_rad',{});

curr = struct('caseid',"",'mach',NaN,'altitude_ft',NaN,'gamma_deg',NaN, ...
    'thrust_lbf_per_engine',NaN,'phase',"unknown",'altitude_m',NaN,'gamma_rad',NaN);
in_case = false;

for i = 1:numel(lines)
    s = strtrim(lines{i});
    if startsWith(s,'CASEID','IgnoreCase',true)
        if in_case && strlength(curr.caseid) > 0
            case_defs(end+1) = finalize_case(curr, ft_to_m); %#ok<AGROW>
        end
        in_case = true;
        curr = struct('caseid',string(strtrim(s)),'mach',NaN,'altitude_ft',NaN, ...
            'gamma_deg',NaN,'thrust_lbf_per_engine',NaN,'phase',infer_phase(s), ...
            'altitude_m',NaN,'gamma_rad',NaN);
        continue;
    end
    if ~in_case, continue; end
    if contains(s,'MACH(1)','IgnoreCase',true) && isnan(curr.mach)
        tok = regexp(s,'MACH\s*\(\s*1\s*\)\s*=\s*([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)','tokens','once');
        if ~isempty(tok), curr.mach = str2double(strrep(tok{1},'D','E')); end
    end
    if contains(s,'ALT(1)','IgnoreCase',true) && isnan(curr.altitude_ft)
        tok = regexp(s,'ALT\s*\(\s*1\s*\)\s*=\s*([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)','tokens','once');
        if ~isempty(tok), curr.altitude_ft = str2double(strrep(tok{1},'D','E')); end
    end
    if contains(s,'GAMMA','IgnoreCase',true) && isnan(curr.gamma_deg)
        tok = regexp(s,'GAMMA\s*=\s*([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)','tokens','once');
        if ~isempty(tok), curr.gamma_deg = str2double(strrep(tok{1},'D','E')); end
    end
    if contains(s,'THSTCJ','IgnoreCase',true) && isnan(curr.thrust_lbf_per_engine)
        tok = regexp(s,'THSTCJ\s*=\s*([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)','tokens','once');
        if ~isempty(tok), curr.thrust_lbf_per_engine = str2double(strrep(tok{1},'D','E')); end
    end
end
if in_case && strlength(curr.caseid) > 0
    case_defs(end+1) = finalize_case(curr, ft_to_m);
end
end

function c = finalize_case(c, ft_to_m)
if isfinite(c.altitude_ft), c.altitude_m = c.altitude_ft*ft_to_m; else, c.altitude_m = NaN; end
if isfinite(c.gamma_deg),   c.gamma_rad  = deg2rad(c.gamma_deg);  else, c.gamma_rad  = NaN; end
end

function phase = infer_phase(s)
u = upper(char(s));
if     contains(u,'TAKEOFF'),  phase = "takeoff";
elseif contains(u,'CLIMB'),    phase = "climb";
elseif contains(u,'CRUISE'),   phase = "cruise";
elseif contains(u,'DESCENT'),  phase = "descent";
elseif contains(u,'APPROACH'), phase = "approach";
else,                           phase = "unknown";
end
end

% ── Table extraction ────────────────────────────────────────────────────────────

function tables = extract_aero_tables(lines)
tables = struct('line',{},'mach',{},'altitude_ft',{},'gamma_deg',{}, ...
    'alpha_deg',{},'CD',{},'CL',{},'CM',{},'config_label',{},'config_priority',{});

current_mach   = NaN;
current_alt    = NaN;
current_gamma  = NaN;
current_config = "";
tcount = 0;

for i = 1:numel(lines)
    s = strtrim(lines{i});

    % New Mach block: read flight conditions only — do NOT reset current_config
    % The config label was already set by the preceding CHARACTERISTICS header
    if contains(s,'MACH') && contains(s,'ALTITUDE') && contains(s,'VELOCITY') && contains(s,'REYNOLDS')
        for j = i+1 : min(i+30, numel(lines))
            sj  = strtrim(lines{j});
            raw = sj;
            if startsWith(raw,'0 '), raw = raw(3:end); end
            vals = str2num(raw); %#ok<ST2NM>
            if numel(vals) >= 6
                m_cand = vals(1); a_cand = vals(2);
                if m_cand >= 0.05 && m_cand <= 2.5 && a_cand >= 0 && a_cand < 100000
                    current_mach = m_cand;
                    current_alt  = a_cand;
                    break;
                end
            end
        end
    end

    % Configuration label — read from lines following the CHARACTERISTICS header
    if contains(s,'CHARACTERISTICS AT ANGLE OF ATTACK')
        current_config = "";
        for j = i+1:min(i+8, numel(lines))
            sl = strtrim(lines{j});
            if isempty(sl), continue; end
            su = upper(sl);
            if contains(su,'CONFIGURATION') || contains(su,'WING') || contains(su,'BODY')
                current_config = string(sl);
                break;
            end
            if startsWith(sl,'Boeing','IgnoreCase',true) || startsWith(sl,'0 ')
                break;
            end
        end
    end

    % Data table header line
    if contains(s,'ALPHA') && contains(s,'CD') && contains(s,'CL') && contains(s,'CM')
        lbl = upper(char(current_config));
        if contains(lbl,'JET POWER') || contains(lbl,'DYNAMIC'), continue; end

        config_priority = 0;
        if     contains(lbl,'WING-BODY-VERTICAL TAIL-HORIZONTAL TAIL'), config_priority = 5;
        elseif contains(lbl,'WING-BODY-HORIZONTAL TAIL'),                config_priority = 4;
        elseif contains(lbl,'WING-BODY-VERTICAL TAIL'),                  config_priority = 3;
        elseif contains(lbl,'WING-BODY'),                                config_priority = 2;
        elseif contains(lbl,'WING'),                                     config_priority = 1;
        end
        if config_priority < 1, continue; end

        alpha = []; CD = []; CL = []; CM = [];
        k = i + 1;
        while k <= numel(lines)
            dl = strtrim(lines{k});
            if isempty(dl), k = k+1; continue; end
            if strcmp(dl,'0') || strcmp(dl,'1'), k = k+1; continue; end
            if contains(dl,'AUTOMATED STABILITY') || startsWith(dl,'0***') || ...
               contains(dl,'DYNAMIC DERIVATIVES') || ...
               contains(dl,'WING SECTION')         || contains(dl,'HORIZONTAL TAIL SECTION') || ...
               contains(dl,'VERTICAL TAIL SECTION')|| startsWith(dl,'CASEID','IgnoreCase',true) || ...
               contains(dl,'CHARACTERISTICS AT ANGLE OF ATTACK')
                break;
            end
            vals = str2num(dl); %#ok<ST2NM>
            if ~isempty(vals) && numel(vals) >= 4
                a = vals(1); cd = vals(2); cl = vals(3); cm = vals(4);
                if isfinite(a) && abs(a) < 40 && isfinite(cl) && isfinite(cd) && isfinite(cm)
                    alpha(end+1,1) = a; CD(end+1,1) = cd; %#ok<AGROW>
                    CL(end+1,1)   = cl; CM(end+1,1) = cm; %#ok<AGROW>
                end
            end
            k = k + 1;
        end

        if numel(alpha) >= 8
            [alpha_u, ia_u] = unique(alpha, 'stable');
            if numel(alpha_u) < 4, continue; end
            if ~isfinite(current_mach) || current_mach <= 0, continue; end
            tcount = tcount + 1;
            tables(tcount).line            = i;
            tables(tcount).mach            = current_mach;
            tables(tcount).altitude_ft     = current_alt;
            tables(tcount).gamma_deg       = current_gamma;
            tables(tcount).alpha_deg       = alpha_u;
            tables(tcount).CD              = CD(ia_u);
            tables(tcount).CL              = CL(ia_u);
            tables(tcount).CM              = CM(ia_u);
            tables(tcount).config_label    = current_config;
            tables(tcount).config_priority = config_priority;
        end
    end
end

fprintf('Raw tables found: %d\n', tcount);
tables = deduplicate_tables(tables);
end

function tables_out = deduplicate_tables(tables_in)
if isempty(tables_in), tables_out = tables_in; return; end
mach_vals   = round([tables_in.mach], 3);
unique_mach = unique(mach_vals);
keep = false(numel(tables_in),1);
for im = 1:numel(unique_mach)
    idx = find(abs(mach_vals - unique_mach(im)) < 0.001);
    [~, best] = max([tables_in(idx).config_priority]);
    keep(idx(best)) = true;
end
tables_out = tables_in(keep);
end