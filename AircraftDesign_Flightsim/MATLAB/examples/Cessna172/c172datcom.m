function [lookup, data] = c172datcom(datcom_file_path,debug_print,model_options)
% C172DATCOM  Build a coefficient lookup from Digital DATCOM .out files.
%
% Usage:
%   lookup = c172datcom(datcom_file_path)
%   lookup = c172datcom(["case010.out","case012.out","case014.out"])
%   [lookup,data] = c172datcom(datcom_file_path,true)
%
% State convention:
%   x = [X; Y; Z; u; v; w; phi; theta; psi; p; q; r]
%
% Controls are resolved by registered name. The legacy first-three ordering
% remains a fallback when no Aircraft object is supplied.

if nargin < 1 || isempty(datcom_file_path)
    error('c172datcom:MissingFile','DATCOM file path required.');
end

if nargin < 2 || isempty(debug_print)
    debug_print = true;
end

if nargin < 3 || isempty(model_options)
    model_options = struct();
end
if ~isstruct(model_options) || ~isscalar(model_options)
    error('c172datcom:InvalidModelOptions', 'model_options must be a scalar struct.');
end

source_files = normalize_datcom_files(datcom_file_path);
data = parse_C172datcom(source_files,model_options);

if debug_print
    print_C172datcom_summary(data,source_files);
end

lookup = @(x,u,geom,varargin) eval_C172datcom(x,u,geom,get_optional_aircraft(varargin),data);

end

% =========================================================================
% Coefficient evaluation
% =========================================================================

function coeff = eval_C172datcom(x,u,geom,aircraft,d)

u_b = double(x(4));
v_b = double(x(5));
w_b = double(x(6));

p_b = double(x(10));
q_b = double(x(11));
r_b = double(x(12));

V = sqrt(u_b^2 + v_b^2 + w_b^2);
V = max(V,1e-9);

altitude_m = max(-double(x(3)),0);
if isa(aircraft,'Aircraft') && isvalid(aircraft)
    [~,speed_of_sound_mps,~,~] = aircraft.get_atmosphere(altitude_m);
else
    [~,speed_of_sound_mps,~,~] = FlightEnvironment.standard_atmosphere(altitude_m);
end
Mach = V/max(speed_of_sound_mps,1e-9);

alpha = atan2(w_b,max(abs(u_b),1e-9));
beta  = asin(max(-1,min(1,v_b/V)));

u_d = double(u(:));

da = get_named_control(u_d,aircraft,'aileron',1);
de = get_named_control(u_d,aircraft,'elevator',2);
dr = get_named_control(u_d,aircraft,'rudder',3);

[~,bref,cref] = get_reference_geometry_safe(geom);

p_hat = p_b * bref / (2*V);
q_hat = q_b * cref / (2*V);
r_hat = r_b * bref / (2*V);

CL_base = interpolate_database_ragged_3d(d.mach,d.altitude_m,d.alpha_grid,d.CL_curve_grid, Mach,altitude_m,alpha);
CD_base = interpolate_database_ragged_3d(d.mach,d.altitude_m,d.alpha_grid,d.CD_curve_grid, Mach,altitude_m,alpha);
CM_base = interpolate_database_ragged_3d(d.mach,d.altitude_m,d.alpha_grid,d.CM_curve_grid, Mach,altitude_m,alpha);

CLA = interpolate_database_2d(d.mach,d.altitude_m,d.CLA_grid,Mach,altitude_m);
CMA = interpolate_database_2d(d.mach,d.altitude_m,d.CMA_grid,Mach,altitude_m);
CYB = interpolate_database_2d(d.mach,d.altitude_m,d.CYB_grid,Mach,altitude_m);
CNB = interpolate_database_2d(d.mach,d.altitude_m,d.CNB_grid,Mach,altitude_m);
CLB = interpolate_database_2d(d.mach,d.altitude_m,d.CLB_grid,Mach,altitude_m);

CLQ = interpolate_database_2d(d.mach,d.altitude_m,d.CLQ_grid,Mach,altitude_m);
CMQ = interpolate_database_2d(d.mach,d.altitude_m,d.CMQ_grid,Mach,altitude_m);
CLP = interpolate_database_2d(d.mach,d.altitude_m,d.CLP_grid,Mach,altitude_m);
CYP = interpolate_database_2d(d.mach,d.altitude_m,d.CYP_grid,Mach,altitude_m);
CNP = interpolate_database_2d(d.mach,d.altitude_m,d.CNP_grid,Mach,altitude_m);
CNR = interpolate_database_2d(d.mach,d.altitude_m,d.CNR_grid,Mach,altitude_m);
CLR = interpolate_database_2d(d.mach,d.altitude_m,d.CLR_grid,Mach,altitude_m);

coeff = struct();

% Longitudinal
coeff.CL = CL_base + d.CLDE * de + CLQ  * q_hat;

coeff.CD = CD_base + d.CD0_add + d.CD_alpha_extra * alpha^2;
coeff.CD_datcom = CD_base;
coeff.CD0_add = d.CD0_add;
coeff.CD_alpha_extra = d.CD_alpha_extra;

coeff.Cm = CM_base + d.CMDE * de + CMQ  * q_hat;

% Lateral-directional
coeff.CY = CYB  * beta + d.CYDR * dr + CYP  * p_hat;

coeff.Cl = CLB  * beta + d.CLDA * da + CLP  * p_hat + CLR  * r_hat;

coeff.Cn = CNB  * beta + d.CNDA * da + d.CNDR * dr + CNP  * p_hat + CNR  * r_hat;

coeff.CL = max(min(coeff.CL,d.CL_max_abs),-d.CL_max_abs);

% Useful diagnostic fields
coeff.alpha_rad = alpha;
coeff.beta_rad  = beta;
coeff.V = V;
coeff.Mach = Mach;
coeff.altitude_m = altitude_m;
coeff.database_Mach = clamp_to_range(Mach,d.mach);
coeff.database_altitude_m = clamp_to_range(altitude_m,d.altitude_m);
alpha_min_local = interpolate_database_2d(d.mach,d.altitude_m,d.alpha_min_grid,Mach,altitude_m);
alpha_max_local = interpolate_database_2d(d.mach,d.altitude_m,d.alpha_max_grid,Mach,altitude_m);
coeff.database_alpha_min_rad = alpha_min_local;
coeff.database_alpha_max_rad = alpha_max_local;
coeff.database_alpha_rad = max(alpha_min_local,min(alpha,alpha_max_local));
coeff.database_alpha_clamped = abs(coeff.database_alpha_rad-alpha) > 1e-12;
coeff.p_hat = p_hat;
coeff.q_hat = q_hat;
coeff.r_hat = r_hat;

% Derivative metadata
coeff.CLA = CLA;
coeff.CMA = CMA;
coeff.CYB = CYB;
coeff.CNB = CNB;
coeff.CLB = CLB;

coeff.CLQ = CLQ;
coeff.CMQ = CMQ;
coeff.CLP = CLP;
coeff.CYP = CYP;
coeff.CNP = CNP;
coeff.CNR = CNR;
coeff.CLR = CLR;

coeff.CLDE = d.CLDE;
coeff.CMDE = d.CMDE;
coeff.CLDA = d.CLDA;
coeff.CNDA = d.CNDA;
coeff.CYDR = d.CYDR;
coeff.CNDR = d.CNDR;

% Conventional aliases
coeff.Clb = CLB;
coeff.Cnb = CNB;
coeff.Cyb = CYB;

coeff.CLq = CLQ;
coeff.Cmq = CMQ;

coeff.Clp = CLP;
coeff.Cyp = CYP;
coeff.Cnp = CNP;
coeff.Cnr = CNR;
coeff.Clr = CLR;

coeff.takeoff_params = d.takeoff_params;
coeff.landing_params = d.landing_params;
coeff.control_effects_included = true;
coeff.output_axes = "body";
coeff.moment_reference_point_body_m = geom.get_reference_point();

end

function value = get_named_control(u,aircraft,name,legacy_index)
value = 0;
if isa(aircraft,'Aircraft') && isvalid(aircraft)
    index = aircraft.control.get_index_by_name(name);
    if index > 0 && index <= numel(u)
        value = u(index);
    end
elseif legacy_index <= numel(u)
    value = u(legacy_index);
end
end

function aircraft = get_optional_aircraft(arguments)
aircraft = [];
if ~isempty(arguments)
    aircraft = arguments{1};
end
end

% =========================================================================
% DATCOM parser
% =========================================================================

function data = parse_C172datcom(source_files,model_options)

cases = repmat(empty_datcom_case(),0,1);

for file_index = 1:numel(source_files)
    file_cases = parse_datcom_file(source_files{file_index});
    cases = [cases; file_cases(:)]; %#ok<AGROW>
end

if isempty(cases)
    error('c172datcom:NoStaticTable', 'No full-configuration static CL/CD/CM tables were found.');
end

mach = unique([cases.mach].','sorted');
altitude_m = unique([cases.altitude_m].','sorted');
original_alpha_counts = arrayfun(@(item) numel(item.alpha),cases).';
alpha = unique(vertcat(cases.alpha),'sorted');
if any(original_alpha_counts < 5)
    error('c172datcom:InsufficientAlphaGrid', 'Every DATCOM condition must contain at least five alpha values.');
end

n_mach = numel(mach);
n_altitude = numel(altitude_m);
n_alpha = numel(alpha);

CL_grid = nan(n_mach,n_altitude,n_alpha);
CD_grid = nan(n_mach,n_altitude,n_alpha);
CM_grid = nan(n_mach,n_altitude,n_alpha);
alpha_grid = cell(n_mach,n_altitude);
CL_curve_grid = cell(n_mach,n_altitude);
CD_curve_grid = cell(n_mach,n_altitude);
CM_curve_grid = cell(n_mach,n_altitude);
alpha_min_grid = nan(n_mach,n_altitude);
alpha_max_grid = nan(n_mach,n_altitude);

static_names = {'CLA','CMA','CYB','CNB','CLB'};
dynamic_names = {'CLQ','CMQ','CLP','CYP','CNP','CNR','CLR'};

for name_index = 1:numel(static_names)
    static_grids.(static_names{name_index}) = nan(n_mach,n_altitude);
end

for name_index = 1:numel(dynamic_names)
    dynamic_grids.(dynamic_names{name_index}) = nan(n_mach,n_altitude);
    raw_dynamic_grids.(dynamic_names{name_index}) = nan(n_mach,n_altitude);
end

dynamic_tables = cell(n_mach,n_altitude);
source_grid = strings(n_mach,n_altitude);

for case_index = 1:numel(cases)
    i_mach = find(abs(mach-cases(case_index).mach) < 1e-10,1);
    i_altitude = find(abs(altitude_m-cases(case_index).altitude_m) < 1e-7,1);

    if ~isempty(alpha_grid{i_mach,i_altitude})
        error('c172datcom:DuplicateCondition', 'Duplicate DATCOM condition at Mach %.4f and altitude %.3f m.', cases(case_index).mach,cases(case_index).altitude_m);
    end

    case_alpha = cases(case_index).alpha(:);
    case_CL = cases(case_index).CL(:);
    case_CD = cases(case_index).CD(:);
    case_CM = cases(case_index).CM(:);
    alpha_grid{i_mach,i_altitude} = case_alpha;
    CL_curve_grid{i_mach,i_altitude} = case_CL;
    CD_curve_grid{i_mach,i_altitude} = case_CD;
    CM_curve_grid{i_mach,i_altitude} = case_CM;
    alpha_min_grid(i_mach,i_altitude) = min(case_alpha);
    alpha_max_grid(i_mach,i_altitude) = max(case_alpha);

    for i_alpha = 1:numel(case_alpha)
        [difference,union_index] = min(abs(alpha-case_alpha(i_alpha)));
        if difference > 1e-10
            error('c172datcom:AlphaUnionFailure', 'Could not map a case onto the alpha union.');
        end
        CL_grid(i_mach,i_altitude,union_index) = case_CL(i_alpha);
        CD_grid(i_mach,i_altitude,union_index) = case_CD(i_alpha);
        CM_grid(i_mach,i_altitude,union_index) = case_CM(i_alpha);
    end

    for name_index = 1:numel(static_names)
        name = static_names{name_index};
        static_grids.(name)(i_mach,i_altitude) = cases(case_index).(name);
    end

    for name_index = 1:numel(dynamic_names)
        name = dynamic_names{name_index};
        dynamic_grids.(name)(i_mach,i_altitude) = cases(case_index).(name);
        raw_dynamic_grids.(name)(i_mach,i_altitude) = cases(case_index).raw_dynamic_derivatives.(name);
    end

    dynamic_tables{i_mach,i_altitude} = cases(case_index).dynamic_table;
    source_grid(i_mach,i_altitude) = string(cases(case_index).source_file);
end

if any(cellfun(@isempty,alpha_grid(:)))
    error('c172datcom:IncompleteConditionGrid', 'The DATCOM files do not form a complete Mach-altitude grid.');
end

reference = vertcat(cases.reference_values);
reference0 = reference(1,:);
if any(max(abs(reference-reference0),[],1) > 1e-6)
    error('c172datcom:InconsistentReferenceGeometry', 'DATCOM reference area, lengths, or moment center differ between files.');
end

data = struct();
data.source_files = string(source_files(:));
data.source_grid = source_grid;
data.mach = mach;
data.altitude_m = altitude_m;
data.alpha = alpha;
data.CL_grid = CL_grid;
data.CD_grid = CD_grid;
data.CM_grid = CM_grid;
data.alpha_grid = alpha_grid;
data.CL_curve_grid = CL_curve_grid;
data.CD_curve_grid = CD_curve_grid;
data.CM_curve_grid = CM_curve_grid;
data.alpha_min_grid = alpha_min_grid;
data.alpha_max_grid = alpha_max_grid;
data.CL = CL_curve_grid{1,1};
data.CD = CD_curve_grid{1,1};
data.CM = CM_curve_grid{1,1};
data.database_size = [n_mach,n_altitude,n_alpha];
data.original_alpha_counts = original_alpha_counts;
data.ragged_alpha_grid = any(original_alpha_counts ~= n_alpha);
data.alpha_intersection_applied = false;

for name_index = 1:numel(static_names)
    name = static_names{name_index};
    data.([name '_grid']) = static_grids.(name);
    data.(name) = static_grids.(name)(1,1);
end

for name_index = 1:numel(dynamic_names)
    name = dynamic_names{name_index};
    data.([name '_grid']) = dynamic_grids.(name);
    data.(name) = dynamic_grids.(name)(1,1);
end

data.raw_dynamic_derivatives = raw_dynamic_grids;
data.dynamic_table = dynamic_tables;
data.reference_area_m2 = reference0(1);
data.reference_chord_m = reference0(2);
data.reference_span_m = reference0(3);
data.moment_reference_horizontal_m = reference0(4);
data.moment_reference_vertical_m = reference0(5);

% Manual control derivatives
data.CLDE =  0.30;
data.CMDE = -1.60;

data.CLDA =  0.075;
data.CNDA =  0.008;

data.CYDR =  0.20;
data.CNDR = -0.075;

% Extra drag / coefficient limits
data.CD0_add = 0.0;
data.CD_alpha_extra = 0.0;
data.CL_max_abs = 1.9;

data = apply_model_options(data,model_options);

% Takeoff / landing model parameters
data.takeoff_params = struct( ...
    'mu_rolling',0.03, ...
    'mu_braking',0.35, ...
    'safety_factor',1.15, ...
    'CLmax_takeoff',1.6, ...
    'VR_to_Vs_ratio',1.10, ...
    'V2_to_Vs_ratio',1.20, ...
    'screen_height_takeoff_ft',50, ...
    'rotation_alpha_deg',8, ...
    'reaction_time_s',1.0);

data.landing_params = struct( ...
    'mu_braking',0.50, ...
    'approach_angle_deg',3.0, ...
    'safety_factor',1.67, ...
    'CLmax_landing',1.9, ...
    'Vapp_to_Vs_ratio',1.30, ...
    'screen_height_landing_ft',50, ...
    'flare_height_m',3.0, ...
    'idle_throttle',0.05, ...
    'use_idle_thrust',true);

end

function data = apply_model_options(data,options)

allowed_scalar_fields = {'CLDE','CMDE','CLDA','CNDA','CYDR','CNDR', 'CD0_add','CD_alpha_extra','CL_max_abs'};
option_names = fieldnames(options);

for k = 1:numel(option_names)
    name = option_names{k};
    if ~ismember(name,allowed_scalar_fields)
        error('c172datcom:UnknownModelOption', 'Unknown model option "%s".',name);
    end
    value = options.(name);
    if ~isscalar(value) || ~isreal(value) || ~isfinite(value)
        error('c172datcom:InvalidModelOption', 'Model option "%s" must be a finite real scalar.',name);
    end
    data.(name) = value;
end

if data.CL_max_abs <= 0
    error('c172datcom:InvalidCLLimit', 'CL_max_abs must be positive.');
end

data.model_options = options;

end

function file_cases = parse_datcom_file(filepath)

txt = fileread(filepath);
lines = regexp(txt,'\r\n|\n|\r','split');
static_headers = [];

for line_index = 1:numel(lines)
    if is_static_table_header(lines{line_index})
        context = strjoin(lines(max(1,line_index-25):line_index),newline);
        if contains(context, 'WING-BODY-VERTICAL TAIL-HORIZONTAL TAIL CONFIGURATION')
            static_headers(end+1,1) = line_index; %#ok<AGROW>
        end
    end
end

if isempty(static_headers)
    error('c172datcom:NoStaticTable', 'No full-configuration static table found in %s.',filepath);
end

file_cases = repmat(empty_datcom_case(),numel(static_headers),1);

for table_index = 1:numel(static_headers)
    header_index = static_headers(table_index);
    if table_index < numel(static_headers)
        segment_end = static_headers(table_index+1)-1;
    else
        segment_end = numel(lines);
    end
    segment = lines(max(1,header_index-25):segment_end);

    [alpha_deg,CL,CD,CM,CLA,CMA,CYB,CNB,CLB] = parse_static_table(segment);
    [CLQ,CMQ,CLP,CYP,CNP,CNR,CLR,dyn_table] = parse_dynamic_table(segment);
    condition = parse_flight_condition(lines,header_index,filepath);

    [alpha_unique,unique_index] = unique(alpha_deg(:),'sorted');
    case_data = empty_datcom_case();
    case_data.source_file = filepath;
    case_data.mach = condition.mach;
    case_data.altitude_m = condition.altitude_m;
    case_data.reference_values = condition.reference_values;
    case_data.alpha = deg2rad(alpha_unique);
    case_data.CL = CL(unique_index);
    case_data.CD = CD(unique_index);
    case_data.CM = CM(unique_index);

    if isnan(CLA), CLA = estimate_slope(case_data.alpha,case_data.CL,0); end
    if isnan(CMA), CMA = estimate_slope(case_data.alpha,case_data.CM,0); end
    if isnan(CYB), CYB = -0.25; end
    if isnan(CNB), CNB =  0.08; end
    if isnan(CLB), CLB = -0.18; end

    case_data.CLA = CLA;
    case_data.CMA = CMA;
    case_data.CYB = CYB;
    case_data.CNB = CNB;
    case_data.CLB = CLB;

    raw = struct('CLQ',CLQ,'CMQ',CMQ,'CLP',CLP,'CYP',CYP, 'CNP',CNP,'CNR',CNR,'CLR',CLR);
    case_data.raw_dynamic_derivatives = raw;
    case_data.CLQ = sanitize_dynamic_derivative(CLQ,'CLQ');
    case_data.CMQ = sanitize_dynamic_derivative(CMQ,'CMQ');
    case_data.CLP = sanitize_dynamic_derivative(CLP,'CLP');
    case_data.CYP = sanitize_dynamic_derivative(CYP,'CYP');
    case_data.CNP = sanitize_dynamic_derivative(CNP,'CNP');
    case_data.CNR = sanitize_dynamic_derivative(CNR,'CNR');
    case_data.CLR = sanitize_dynamic_derivative(CLR,'CLR');
    case_data.dynamic_table = dyn_table;
    file_cases(table_index) = case_data;
end

end

function case_data = empty_datcom_case()

case_data = struct( ...
    'source_file','', ...
    'mach',NaN, ...
    'altitude_m',NaN, ...
    'reference_values',nan(1,5), ...
    'alpha',[], ...
    'CL',[], ...
    'CD',[], ...
    'CM',[], ...
    'CLA',NaN, ...
    'CMA',NaN, ...
    'CYB',NaN, ...
    'CNB',NaN, ...
    'CLB',NaN, ...
    'CLQ',NaN, ...
    'CMQ',NaN, ...
    'CLP',NaN, ...
    'CYP',NaN, ...
    'CNP',NaN, ...
    'CNR',NaN, ...
    'CLR',NaN, ...
    'raw_dynamic_derivatives',struct(), ...
    'dynamic_table',struct());

end

function condition = parse_flight_condition(lines,header_index,filepath)

for line_index = header_index-1:-1:max(1,header_index-20)
    values = read_numeric_values(lines{line_index});
    if numel(values) >= 12 && abs(values(1)) < 1e-12 && values(2) >= 0 && values(2) < 5 && values(8) > 0
        condition = struct();
        condition.mach = values(2);
        condition.altitude_m = values(3);
        condition.reference_values = values(8:12);
        return;
    end
end

error('c172datcom:MissingFlightCondition', 'Could not read flight-condition metadata near line %d of %s.', header_index,filepath);

end

function value = sanitize_dynamic_derivative(value,name)

fallbacks = struct('CLQ',7.35, 'CMQ',-14.8, 'CLP',-0.4752, 'CYP',-0.04625, 'CNP',-0.01073, 'CNR',-0.04790, 'CLR',0.03454);

invalid = ~isfinite(value);
switch name
    case 'CLQ'
        invalid = invalid || abs(value) < 0.1 || abs(value) > 30;
    case 'CMQ'
        invalid = invalid || value >= 0 || abs(value) > 100;
    case 'CLP'
        invalid = invalid || value >= 0 || abs(value) > 5;
    case 'CYP'
        invalid = invalid || abs(value) > 5;
    case 'CNP'
        invalid = invalid || abs(value) > 2;
    case 'CNR'
        invalid = invalid || value >= 0 || abs(value) > 5;
    case 'CLR'
        invalid = invalid || abs(value) > 2;
end

if invalid
    value = fallbacks.(name);
end

end

function tf = is_static_table_header(header)

tf = contains(header,'ALPHA') && contains(header,'CD') && contains(header,'CL') && contains(header,'CM') && contains(header,'XCP') && ~contains(header,'CLQ') && ~contains(header,'Q/QINF');

end

% =========================================================================
% Static table parsing
% =========================================================================

function [best_alpha,best_CL,best_CD,best_CM,best_CLA,best_CMA,best_CYB,best_CNB,best_CLB] = parse_static_table(lines)

best_alpha = [];
best_CL = [];
best_CD = [];
best_CM = [];

best_CLA = NaN;
best_CMA = NaN;
best_CYB = NaN;
best_CNB = NaN;
best_CLB = NaN;

best_score = -inf;

for i = 1:numel(lines)

    header = lines{i};

    is_static_header = contains(header,'ALPHA') && contains(header,'CD') && contains(header,'CL') && contains(header,'CM') && contains(header,'XCP') && ~contains(header,'CLQ') && ~contains(header,'Q/QINF');

    if ~is_static_header
        continue;
    end

    % Prefer full aircraft configuration tables.
    context0 = max(1,i-25);
    context = strjoin(lines(context0:i),newline);

    full_config_bonus = 0;

    if contains(context,'WING-BODY-VERTICAL TAIL-HORIZONTAL TAIL CONFIGURATION')
        full_config_bonus = 1000;
    end

    alpha = [];
    CL = [];
    CD = [];
    CM = [];

    CLA_vals = [];
    CMA_vals = [];
    CYB_vals = [];
    CNB_vals = [];
    CLB_vals = [];

    for k = i+1:min(i+80,numel(lines))

        line = strtrim(lines{k});

        if isempty(line)
            continue;
        end

        if contains(line,'Q/QINF') || contains(line,'DYNAMIC DERIVATIVES') || contains(line,'AUTOMATED STABILITY') || contains(line,'CONTROL DERIVATIVES') || startsWith(line,'0***')
            break;
        end

        tokens = regexp(line,'\S+','match');

        if numel(tokens) < 4
            continue;
        end

        a  = str2double(tokens{1});
        cd = str2double(tokens{2});
        cl = str2double(tokens{3});
        cm = str2double(tokens{4});

        if isfinite(a) && isfinite(cd) && isfinite(cl) && isfinite(cm) && abs(a) <= 35

            alpha(end+1,1) = a; %#ok<AGROW>
            CD(end+1,1) = cd; %#ok<AGROW>
            CL(end+1,1) = cl; %#ok<AGROW>
            CM(end+1,1) = cm; %#ok<AGROW>

            if numel(tokens) >= 8
                val = str2double(tokens{8});
                if isfinite(val), CLA_vals(end+1,1) = val; end %#ok<AGROW>
            end

            if numel(tokens) >= 9
                val = str2double(tokens{9});
                if isfinite(val), CMA_vals(end+1,1) = val; end %#ok<AGROW>
            end

            % DATCOM often leaves CYB/CNB blank on most rows.
            % Full rows have 12 tokens:
            % alpha CD CL CM CN CA XCP CLA CMA CYB CNB CLB
            % Short rows have 10 tokens:
            % alpha CD CL CM CN CA XCP CLA CMA CLB
            if numel(tokens) >= 12
                val = str2double(tokens{10});
                if isfinite(val), CYB_vals(end+1,1) = val; end %#ok<AGROW>

                val = str2double(tokens{11});
                if isfinite(val), CNB_vals(end+1,1) = val; end %#ok<AGROW>

                val = str2double(tokens{12});
                if isfinite(val), CLB_vals(end+1,1) = val; end %#ok<AGROW>
            elseif numel(tokens) >= 10
                val = str2double(tokens{10});
                if isfinite(val), CLB_vals(end+1,1) = val; end %#ok<AGROW>
            end
        end
    end

    score = full_config_bonus + numel(alpha);

    if numel(alpha) >= 5 && score > best_score

        best_score = score;

        best_alpha = alpha;
        best_CL = CL;
        best_CD = CD;
        best_CM = CM;

        best_CLA = select_value_near_zero(best_alpha,CLA_vals);
        best_CMA = select_value_near_zero(best_alpha,CMA_vals);

        if ~isempty(CYB_vals)
            best_CYB = CYB_vals(1);
        else
            best_CYB = NaN;
        end

        if ~isempty(CNB_vals)
            best_CNB = CNB_vals(1);
        else
            best_CNB = NaN;
        end

        if ~isempty(CLB_vals)
            best_CLB = CLB_vals(1);
        else
            best_CLB = NaN;
        end
    end
end

end

% =========================================================================
% Dynamic table parsing
% =========================================================================

function [CLQ,CMQ,CLP,CYP,CNP,CNR,CLR,dyn_table] = parse_dynamic_table(lines)

CLQ = NaN;
CMQ = NaN;
CLP = NaN;
CYP = NaN;
CNP = NaN;
CNR = NaN;
CLR = NaN;

dyn_table = struct('alpha_deg',[], 'CLQ',[], 'CMQ',[], 'CLAD',[], 'CMAD',[], 'CLP',[], 'CYP',[], 'CNP',[], 'CNR',[], 'CLR',[]);

best_score = -inf;
best_table = [];

for i = 1:numel(lines)

    header = lines{i};

    is_dynamic_header = contains(header,'ALPHA') && contains(header,'CLQ') && contains(header,'CMQ') && contains(header,'CLP') && contains(header,'CNR') && contains(header,'CLR');

    if ~is_dynamic_header
        continue;
    end

    context0 = max(1,i-25);
    context = strjoin(lines(context0:i),newline);

    full_config_bonus = 0;

    if contains(context,'WING-BODY-VERTICAL TAIL-HORIZONTAL TAIL CONFIGURATION')
        full_config_bonus = 1000;
    end

    tab = struct('alpha_deg',[], 'CLQ',[], 'CMQ',[], 'CLAD',[], 'CMAD',[], 'CLP',[], 'CYP',[], 'CNP',[], 'CNR',[], 'CLR',[]);

    current_CLQ = NaN;
    current_CMQ = NaN;

    for k = i+1:min(i+80,numel(lines))

        line = strtrim(lines{k});

        if isempty(line)
            continue;
        end

        if contains(line,'AUTOMATED STABILITY') || contains(line,'CONTROL DERIVATIVES') || contains(line,'CHARACTERISTICS AT ANGLE') || contains(line,'ALPHA') || startsWith(line,'0***')
            break;
        end

        vals = read_numeric_values(line);

        if isempty(vals)
            continue;
        end

        a = vals(1);

        if ~isfinite(a) || abs(a) > 35
            continue;
        end

        if numel(vals) >= 10
            % Full row:
            % alpha CLQ CMQ CLAD CMAD CLP CYP CNP CNR CLR
            clq  = vals(2);
            cmq  = vals(3);
            clad = vals(4);
            cmad = vals(5);
            clp  = vals(6);
            cyp  = vals(7);
            cnp  = vals(8);
            cnr  = vals(9);
            clr  = vals(10);

            if isfinite(clq), current_CLQ = clq; end
            if isfinite(cmq), current_CMQ = cmq; end

        elseif numel(vals) >= 8
            % DATCOM omits repeated CLQ/CMQ in later rows:
            % alpha CLAD CMAD CLP CYP CNP CNR CLR
            clq  = current_CLQ;
            cmq  = current_CMQ;
            clad = vals(2);
            cmad = vals(3);
            clp  = vals(4);
            cyp  = vals(5);
            cnp  = vals(6);
            cnr  = vals(7);
            clr  = vals(8);

        else
            continue;
        end

        tab.alpha_deg(end+1,1) = a; %#ok<AGROW>
        tab.CLQ(end+1,1) = clq; %#ok<AGROW>
        tab.CMQ(end+1,1) = cmq; %#ok<AGROW>
        tab.CLAD(end+1,1) = clad; %#ok<AGROW>
        tab.CMAD(end+1,1) = cmad; %#ok<AGROW>
        tab.CLP(end+1,1) = clp; %#ok<AGROW>
        tab.CYP(end+1,1) = cyp; %#ok<AGROW>
        tab.CNP(end+1,1) = cnp; %#ok<AGROW>
        tab.CNR(end+1,1) = cnr; %#ok<AGROW>
        tab.CLR(end+1,1) = clr; %#ok<AGROW>
    end

    score = full_config_bonus + numel(tab.alpha_deg);

    if numel(tab.alpha_deg) >= 5 && score > best_score
        best_score = score;
        best_table = tab;
    end
end

if isempty(best_table)
    return;
end

dyn_table = best_table;

% Pick derivatives near alpha = 0 deg.
CLQ = select_near_alpha(dyn_table.alpha_deg,dyn_table.CLQ,0);
CMQ = select_near_alpha(dyn_table.alpha_deg,dyn_table.CMQ,0);
CLP = select_near_alpha(dyn_table.alpha_deg,dyn_table.CLP,0);
CYP = select_near_alpha(dyn_table.alpha_deg,dyn_table.CYP,0);
CNP = select_near_alpha(dyn_table.alpha_deg,dyn_table.CNP,0);
CNR = select_near_alpha(dyn_table.alpha_deg,dyn_table.CNR,0);
CLR = select_near_alpha(dyn_table.alpha_deg,dyn_table.CLR,0);

end

% =========================================================================
% Utilities
% =========================================================================

function source_files = normalize_datcom_files(input_value)

if ischar(input_value)
    input_value = string(input_value);
end

if isstring(input_value) && isscalar(input_value) && isfolder(input_value)
    listing = dir(fullfile(char(input_value),'*.out'));
    source_files = arrayfun(@(item) fullfile(item.folder,item.name),listing,'UniformOutput',false);
elseif isstring(input_value)
    source_files = cellstr(input_value(:));
elseif iscell(input_value)
    source_files = cell(size(input_value(:)));
    for file_index = 1:numel(input_value)
        if ~(ischar(input_value{file_index}) || (isstring(input_value{file_index}) && isscalar(input_value{file_index})))
            error('c172datcom:InvalidFileList', 'Each DATCOM file path must be a character vector or string scalar.');
        end
        source_files{file_index} = char(string(input_value{file_index}));
    end
else
    error('c172datcom:InvalidFileList', 'Input must be a file path, folder path, string array, or cell array.');
end

if isempty(source_files)
    error('c172datcom:NoFiles', 'No DATCOM .out files were supplied or found in the folder.');
end

for file_index = 1:numel(source_files)
    source_files{file_index} = strtrim(source_files{file_index});
    if ~exist(source_files{file_index},'file')
        error('c172datcom:FileNotFound', 'File not found: %s',source_files{file_index});
    end
end

[~,order] = sort(lower(string(source_files)));
source_files = source_files(order);

end

function value = interpolate_database_ragged_3d(mach_axis,altitude_axis,alpha_grid,curve_grid, mach_query,altitude_query,alpha_query)

mach_query = clamp_to_range(mach_query,mach_axis);
altitude_query = clamp_to_range(altitude_query,altitude_axis);
value_by_condition = zeros(numel(mach_axis),numel(altitude_axis));

for i_mach = 1:numel(mach_axis)
    for i_altitude = 1:numel(altitude_axis)
        local_alpha = alpha_grid{i_mach,i_altitude}(:);
        local_curve = curve_grid{i_mach,i_altitude}(:);
        local_query = clamp_to_range(alpha_query,local_alpha);
        value_by_condition(i_mach,i_altitude) = interpolate_vector(local_alpha,local_curve,local_query);
    end
end

value = interpolate_database_2d(mach_axis,altitude_axis,value_by_condition,mach_query,altitude_query);

end

function value = interpolate_database_2d(mach_axis,altitude_axis,grid,mach_query,altitude_query)

mach_query = clamp_to_range(mach_query,mach_axis);
altitude_query = clamp_to_range(altitude_query,altitude_axis);

value_by_mach = zeros(numel(mach_axis),1);
for i_mach = 1:numel(mach_axis)
    value_by_mach(i_mach) = interpolate_vector(altitude_axis,reshape(grid(i_mach,:),[],1),altitude_query);
end

value = interpolate_vector(mach_axis,value_by_mach,mach_query);

end

function value = interpolate_vector(axis,values,query)

axis = axis(:);
values = values(:);

if numel(axis) == 1
    value = values(1);
else
    value = interp1(axis,values,query,'linear');
end

end

function value = clamp_to_range(value,axis)

axis = axis(:);
value = max(min(value,max(axis)),min(axis));

end

function vals = read_numeric_values(line)

line = strrep(line,'D','E');
line = strrep(line,'d','e');

numstr = regexp(line,'[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?','match');

vals = zeros(1,numel(numstr));

for i = 1:numel(numstr)
    vals(i) = str2double(numstr{i});
end

end

function value = select_near_alpha(alpha_deg,values,alpha_target)

value = NaN;

if isempty(alpha_deg) || isempty(values)
    return;
end

alpha_deg = alpha_deg(:);
values = values(:);

valid = isfinite(alpha_deg) & isfinite(values);

if ~any(valid)
    return;
end

a = alpha_deg(valid);
v = values(valid);

[~,idx] = min(abs(a-alpha_target));
value = v(idx);

end

function value = select_value_near_zero(alpha_deg,values)

% This helper is intentionally conservative. For DATCOM static table
% derivatives, some fields are only printed on sparse rows. Use the first
% finite value if the lengths do not match the alpha table.

value = NaN;

if isempty(values)
    return;
end

values = values(:);
values = values(isfinite(values));

if isempty(values)
    return;
end

if numel(values) == numel(alpha_deg)
    [~,idx] = min(abs(alpha_deg(:)));
    value = values(idx);
else
    value = values(1);
end

end

function slope = estimate_slope(x,y,x0)

x = x(:);
y = y(:);

valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

if numel(x) < 2
    slope = NaN;
    return;
end

[~,idx] = min(abs(x-x0));

if idx == 1
    idx1 = 1;
    idx2 = 2;
elseif idx == numel(x)
    idx1 = numel(x)-1;
    idx2 = numel(x);
else
    idx1 = idx-1;
    idx2 = idx+1;
end

dx = x(idx2)-x(idx1);

if abs(dx) < 1e-12
    slope = NaN;
else
    slope = (y(idx2)-y(idx1))/dx;
end

end

function [Sref,bref,cref] = get_reference_geometry_safe(geom)

Sref = 1;
bref = 1;
cref = 1;

try
    if ismethod(geom,'get_reference_geometry')
        [Sref,bref,cref] = geom.get_reference_geometry();

        if Sref > 0 && bref > 0 && cref > 0
            return;
        end
    end
catch
end

Sref = get_geom_value(geom,{'ref_area','wing_area','S_ref','S'},1);
bref = get_geom_value(geom,{'ref_span','wing_span','b_ref','b'},1);
cref = get_geom_value(geom,{'ref_chord','mean_aerodynamic_chord','c_ref','c'},1);

if cref <= 0 && Sref > 0 && bref > 0
    cref = Sref / bref;
end

Sref = max(Sref,1e-9);
bref = max(bref,1e-9);
cref = max(cref,1e-9);

end

function value = get_geom_value(geom,names,default_value)

value = default_value;

for i = 1:numel(names)

    name = names{i};

    try
        if isprop(geom,name)
            temp = geom.(name);

            if ~isempty(temp) && isnumeric(temp) && isscalar(temp) && isfinite(temp)
                value = temp;
                return;
            end
        end
    catch
    end

    try
        if isstruct(geom) && isfield(geom,name)
            temp = geom.(name);

            if ~isempty(temp) && isnumeric(temp) && isscalar(temp) && isfinite(temp)
                value = temp;
                return;
            end
        end
    catch
    end
end

end

function print_C172datcom_summary(data,source_files)

fprintf('\n=== C172 DATCOM PARSED DATA ===\n');
fprintf('Files             : %d\n',numel(source_files));
for file_index = 1:numel(source_files)
    fprintf('  %s\n',source_files{file_index});
end
fprintf('Database size     : %d Mach x %d altitude x %d alpha\n', data.database_size(1),data.database_size(2),data.database_size(3));
fprintf('Mach range        : %.3f to %.3f\n',min(data.mach),max(data.mach));
fprintf('Altitude range    : %.1f to %.1f m\n', min(data.altitude_m),max(data.altitude_m));
fprintf('Alpha table points: %d\n',numel(data.alpha));
fprintf('Alpha range       : %.2f to %.2f deg\n', rad2deg(min(data.alpha)),rad2deg(max(data.alpha)));
fprintf('Alpha intersection: %d\n',data.alpha_intersection_applied);
fprintf('Ragged alpha grid : %d\n',data.ragged_alpha_grid);
if data.ragged_alpha_grid
    fprintf('Original alpha pts: %d to %d per condition\n', min(data.original_alpha_counts),max(data.original_alpha_counts));
end
fprintf('CL range          : %.4f to %.4f\n', min(data.CL_grid(:),[],'omitnan'),max(data.CL_grid(:),[],'omitnan'));
fprintf('CD range          : %.4f to %.4f\n', min(data.CD_grid(:),[],'omitnan'),max(data.CD_grid(:),[],'omitnan'));
fprintf('CM range          : %.4f to %.4f\n', min(data.CM_grid(:),[],'omitnan'),max(data.CM_grid(:),[],'omitnan'));
fprintf('Configured CD add : %.5f\n',data.CD0_add);
fprintf('Moment center     : horizontal %.3f m, vertical %.3f m\n', data.moment_reference_horizontal_m,data.moment_reference_vertical_m);

fprintf('\n--- Static derivatives ---\n');
fprintf('CLA  = %+.6f\n',data.CLA);
fprintf('CMA  = %+.6f\n',data.CMA);
fprintf('CYB  = %+.6f\n',data.CYB);
fprintf('CNB  = %+.6f\n',data.CNB);
fprintf('CLB  = %+.6f\n',data.CLB);

fprintf('\n--- Dynamic derivatives near alpha = 0 deg ---\n');
fprintf('CLQ  = %+.6f\n',data.CLQ);
fprintf('CMQ  = %+.6f\n',data.CMQ);
fprintf('CLP  = %+.6f\n',data.CLP);
fprintf('CYP  = %+.6f\n',data.CYP);
fprintf('CNP  = %+.6f\n',data.CNP);
fprintf('CNR  = %+.6f\n',data.CNR);
fprintf('CLR  = %+.6f\n',data.CLR);

if isfield(data,'raw_dynamic_derivatives')
    r = data.raw_dynamic_derivatives;

    fprintf('\n--- Raw dynamic derivative ranges before sanity check ---\n');
    fprintf('CLQ_raw = %+.6f to %+.6f\n',min(r.CLQ(:)),max(r.CLQ(:)));
    fprintf('CMQ_raw = %+.6f to %+.6f\n',min(r.CMQ(:)),max(r.CMQ(:)));
    fprintf('CLP_raw = %+.6f to %+.6f\n',min(r.CLP(:)),max(r.CLP(:)));
    fprintf('CYP_raw = %+.6f to %+.6f\n',min(r.CYP(:)),max(r.CYP(:)));
    fprintf('CNP_raw = %+.6f to %+.6f\n',min(r.CNP(:)),max(r.CNP(:)));
    fprintf('CNR_raw = %+.6f to %+.6f\n',min(r.CNR(:)),max(r.CNR(:)));
    fprintf('CLR_raw = %+.6f to %+.6f\n',min(r.CLR(:)),max(r.CLR(:)));
end

fprintf('================================\n');

end
