function [lookup,data] = b737datcom(datcom_files,verbose,control_model)

if nargin < 2 || isempty(verbose)
    verbose = true;
end
if nargin < 3 || isempty(control_model)
    control_model = default_control_model();
else
    control_model = complete_control_model(control_model);
end

files = normalize_file_list(datcom_files);
static_tables = empty_static_tables();
dynamic_tables = empty_dynamic_tables();

for k = 1:numel(files)
    [static_k,dynamic_k] = parse_datcom_file(files(k));
    static_tables = [static_tables,static_k]; %#ok<AGROW>
    dynamic_tables = [dynamic_tables,dynamic_k]; %#ok<AGROW>
end

static_tables = select_best_tables(static_tables);
dynamic_tables = select_best_tables(dynamic_tables);

if isempty(static_tables)
    error('b737datcom:NoStaticData', 'No usable full-aircraft static DATCOM tables were found.');
end

reference = validate_reference_geometry(static_tables);
grid = build_compatibility_grid(static_tables,dynamic_tables);
altitude_scale_m = max(max(grid.altitude_m)-min(grid.altitude_m),1);
interpolants = build_interpolants(static_tables,dynamic_tables,altitude_scale_m);

data = struct();
data.source_files = files;
data.static_tables = static_tables;
data.dynamic_tables = dynamic_tables;
data.reference = reference;
data.grid = grid;
data.interpolants = interpolants;
data.altitude_scale_m = altitude_scale_m;
data.control = control_model;
data.control_model_status = "provisional_user_replaceable";
data.out_of_range_action = "nearest_with_validity_flag";
data.takeoff_params = struct( ...
    'mu_rolling',0.03,'mu_braking',0.35,'safety_factor',1.15, ...
    'CLmax_takeoff',1.8,'VR_to_Vs_ratio',1.10, ...
    'V2_to_Vs_ratio',1.20,'screen_height_takeoff_ft',35, ...
    'rotation_alpha_deg',8,'reaction_time_s',1.0);
data.landing_params = struct( ...
    'mu_braking',0.45,'approach_angle_deg',3.0, ...
    'safety_factor',1.67,'CLmax_landing',2.2, ...
    'Vapp_to_Vs_ratio',1.30,'screen_height_landing_ft',50, ...
    'flare_height_m',3.0,'idle_throttle',0.05, ...
    'use_idle_thrust',true);

lookup = @(x,u,geom) evaluate_lookup(x,u,geom,data);

if verbose
    print_summary(data);
end
end

function files = normalize_file_list(input_value)
if nargin < 1 || isempty(input_value)
    error('b737datcom:FilesRequired','DATCOM output files are required.');
end

if ischar(input_value) || (isstring(input_value) && isscalar(input_value))
    candidate = string(input_value);
    if isfolder(candidate)
        listing = dir(fullfile(candidate,'*.out'));
        if isempty(listing)
            listing = dir(fullfile(candidate,'*.dat'));
        end
        files = strings(numel(listing),1);
        for k = 1:numel(listing)
            files(k) = string(fullfile(listing(k).folder,listing(k).name));
        end
    else
        files = candidate;
    end
elseif iscell(input_value) || isstring(input_value)
    files = string(input_value(:));
else
    error('b737datcom:InvalidFiles', 'Input must be a path, directory, string array, or cell array.');
end

files = files(strlength(files) > 0);
if isempty(files)
    error('b737datcom:NoFiles','No DATCOM output files were found.');
end
for k = 1:numel(files)
    if ~isfile(files(k))
        error('b737datcom:FileNotFound','File not found: %s',files(k));
    end
end
end

function [static_tables,dynamic_tables] = parse_datcom_file(file_name)
text = fileread(file_name);
text = strrep(text,sprintf('\r\n'),sprintf('\n'));
text = strrep(text,sprintf('\r'),sprintf('\n'));
lines = splitlines(string(text));

static_tables = empty_static_tables();
dynamic_tables = empty_dynamic_tables();

for i = 1:numel(lines)
    line = strtrim(lines(i));
    if contains(line,'CHARACTERISTICS AT ANGLE OF ATTACK AND IN SIDESLIP')
        table = parse_static_section(lines,i,file_name);
        if ~isempty(table)
            static_tables(end+1) = table; %#ok<AGROW>
        end
    elseif strcmpi(line,'DYNAMIC DERIVATIVES')
        table = parse_dynamic_section(lines,i,file_name);
        if ~isempty(table)
            dynamic_tables(end+1) = table; %#ok<AGROW>
        end
    end
end
end

function table = parse_static_section(lines,start_index,file_name)
table = [];
[configuration,case_label,priority] = read_configuration(lines,start_index);
if priority < 1
    return;
end
[condition,condition_line] = read_condition(lines,start_index);
if isempty(condition)
    return;
end

header_index = find_header(lines,condition_line, @(s) contains(s,'ALPHA') && contains(s,'CD') && contains(s,'CL') && contains(s,'CM'),45);
if header_index == 0
    return;
end

alpha = [];
CD = [];
CL = [];
CM = [];
CYB = [];
CNB = [];
CLB = [];
last_CYB = NaN;
last_CNB = NaN;

for k = header_index+1:min(header_index+80,numel(lines))
    line = strtrim(lines(k));
    upper_line = upper(line);
    if contains(upper_line,'ALPHA') && contains(upper_line,'Q/QINF')
        break;
    end
    if contains(upper_line,'AUTOMATED STABILITY') || contains(upper_line,'DYNAMIC DERIVATIVES') || contains(upper_line,'CHARACTERISTICS AT ANGLE')
        break;
    end
    values = numeric_tokens(line);
    if numel(values) < 9
        continue;
    end
    angle = values(1);
    if ~isfinite(angle) || abs(angle) > 45
        continue;
    end
    if numel(values) >= 12
        last_CYB = values(10);
        last_CNB = values(11);
        clb = values(12);
    elseif numel(values) >= 10
        clb = values(10);
    else
        clb = NaN;
    end
    alpha(end+1,1) = angle; %#ok<AGROW>
    CD(end+1,1) = values(2); %#ok<AGROW>
    CL(end+1,1) = values(3); %#ok<AGROW>
    CM(end+1,1) = values(4); %#ok<AGROW>
    CYB(end+1,1) = last_CYB; %#ok<AGROW>
    CNB(end+1,1) = last_CNB; %#ok<AGROW>
    CLB(end+1,1) = clb; %#ok<AGROW>
end

valid = isfinite(alpha) & isfinite(CD) & isfinite(CL) & isfinite(CM);
alpha = alpha(valid);
CD = CD(valid);
CL = CL(valid);
CM = CM(valid);
CYB = CYB(valid);
CNB = CNB(valid);
CLB = CLB(valid);
if numel(alpha) < 4
    return;
end

[alpha,unique_index] = unique(alpha,'stable');
CD = CD(unique_index);
CL = CL(unique_index);
CM = CM(unique_index);
CYB = fill_constant_missing(CYB(unique_index));
CNB = fill_constant_missing(CNB(unique_index));
CLB = fill_linear_missing(alpha,CLB(unique_index));

table = struct( ...
    'source_file',string(file_name), ...
    'configuration',configuration, ...
    'case_label',case_label, ...
    'config_priority',priority, ...
    'mach',condition.mach, ...
    'altitude_m',condition.altitude_m, ...
    'alpha_deg',alpha, ...
    'CD',CD,'CL',CL,'CM',CM, ...
    'CYB',CYB,'CNB',CNB,'CLB',CLB, ...
    'sref_m2',condition.sref_m2, ...
    'cbar_m',condition.cbar_m, ...
    'bref_m',condition.bref_m, ...
    'moment_reference_m',condition.moment_reference_m);
end

function table = parse_dynamic_section(lines,start_index,file_name)
table = [];
[configuration,case_label,priority] = read_configuration(lines,start_index);
if priority < 1
    return;
end
[condition,condition_line] = read_condition(lines,start_index);
if isempty(condition)
    return;
end

header_index = find_header(lines,condition_line, @(s) contains(s,'ALPHA') && contains(s,'CLQ') && contains(s,'CMQ') && contains(s,'CNR'),45);
if header_index == 0
    return;
end

alpha = [];
CLQ = [];
CMQ = [];
CLAD = [];
CMAD = [];
CLP = [];
CYP = [];
CNP = [];
CNR = [];
CLR = [];
last_CLQ = NaN;
last_CMQ = NaN;

for k = header_index+1:min(header_index+65,numel(lines))
    line = strtrim(lines(k));
    upper_line = upper(line);
    if contains(upper_line,'WARNING') || contains(upper_line,'AUTOMATED STABILITY') || contains(upper_line,'CONFIGURATION AUXILIARY') || contains(upper_line,'CHARACTERISTICS AT ANGLE')
        break;
    end
    values = numeric_tokens(line);
    if numel(values) >= 10
        angle = values(1);
        last_CLQ = values(2);
        last_CMQ = values(3);
        remainder = values(4:10);
    elseif numel(values) == 8
        angle = values(1);
        remainder = values(2:8);
    else
        continue;
    end
    if ~isfinite(angle) || abs(angle) > 45
        continue;
    end
    alpha(end+1,1) = angle; %#ok<AGROW>
    CLQ(end+1,1) = last_CLQ; %#ok<AGROW>
    CMQ(end+1,1) = last_CMQ; %#ok<AGROW>
    CLAD(end+1,1) = remainder(1); %#ok<AGROW>
    CMAD(end+1,1) = remainder(2); %#ok<AGROW>
    CLP(end+1,1) = remainder(3); %#ok<AGROW>
    CYP(end+1,1) = remainder(4); %#ok<AGROW>
    CNP(end+1,1) = remainder(5); %#ok<AGROW>
    CNR(end+1,1) = remainder(6); %#ok<AGROW>
    CLR(end+1,1) = remainder(7); %#ok<AGROW>
end

if numel(alpha) < 4
    return;
end
[alpha,unique_index] = unique(alpha,'stable');
CLQ = fill_constant_missing(CLQ(unique_index));
CMQ = fill_constant_missing(CMQ(unique_index));
CLAD = CLAD(unique_index);
CMAD = CMAD(unique_index);
CLP = CLP(unique_index);
CYP = CYP(unique_index);
CNP = CNP(unique_index);
CNR = CNR(unique_index);
CLR = CLR(unique_index);

table = struct( ...
    'source_file',string(file_name), ...
    'configuration',configuration, ...
    'case_label',case_label, ...
    'config_priority',priority, ...
    'mach',condition.mach, ...
    'altitude_m',condition.altitude_m, ...
    'alpha_deg',alpha, ...
    'CLQ',CLQ,'CMQ',CMQ,'CLAD',CLAD,'CMAD',CMAD, ...
    'CLP',CLP,'CYP',CYP,'CNP',CNP,'CNR',CNR,'CLR',CLR);
end

function [configuration,case_label,priority] = read_configuration(lines,start_index)
configuration = "";
case_label = "";
for k = start_index+1:min(start_index+8,numel(lines))
    line = strtrim(lines(k));
    if contains(upper(line),'CONFIGURATION') && strlength(configuration) == 0
        configuration = line;
    elseif strlength(line) > 0 && ~contains(upper(line),'FLIGHT CONDITIONS') && ~contains(upper(line),'REFERENCE DIMENSIONS') && ~contains(upper(line),'CONFIGURATION')
        case_label = line;
    end
    if contains(upper(line),'FLIGHT CONDITIONS')
        break;
    end
end

label = upper(configuration);
if contains(label,'WING-BODY-VERTICAL TAIL-HORIZONTAL TAIL') || contains(label,'WING-BODY-HORIZONTAL TAIL-VERTICAL TAIL')
    priority = 5;
elseif contains(label,'WING-BODY-HORIZONTAL TAIL')
    priority = 4;
elseif contains(label,'WING-BODY-VERTICAL TAIL')
    priority = 3;
elseif contains(label,'WING-BODY')
    priority = 2;
elseif contains(label,'WING')
    priority = 1;
else
    priority = 0;
end
end

function [condition,line_index] = read_condition(lines,start_index)
condition = [];
line_index = 0;
metric = true;
for k = start_index+1:min(start_index+35,numel(lines))
    line = upper(strtrim(lines(k)));
    if contains(line,'FT/SEC') || contains(line,'FT**2')
        metric = false;
    elseif contains(line,'M/SEC') || contains(line,'M**2')
        metric = true;
    end
    if startsWith(line,'0 ')
        values = numeric_tokens(extractAfter(strtrim(lines(k)),2));
    else
        values = numeric_tokens(line);
    end
    if numel(values) < 11
        continue;
    end
    if values(1) < 0 || values(1) > 5 || values(2) < -5000
        continue;
    end
    length_factor = 1;
    area_factor = 1;
    if ~metric
        length_factor = 0.3048;
        area_factor = 0.09290304;
    end
    condition = struct( ...
        'mach',values(1), ...
        'altitude_m',values(2)*length_factor, ...
        'sref_m2',values(7)*area_factor, ...
        'cbar_m',values(8)*length_factor, ...
        'bref_m',values(9)*length_factor, ...
        'moment_reference_m',[values(10),0,values(11)]*length_factor);
    line_index = k;
    return;
end
end

function index = find_header(lines,start_index,predicate,max_search)
index = 0;
for k = start_index+1:min(start_index+max_search,numel(lines))
    line = upper(strtrim(lines(k)));
    if predicate(line)
        index = k;
        return;
    end
end
end

function values = numeric_tokens(line)
tokens = regexp(char(line), '[-+]?(?:(?:\d+\.?\d*)|(?:\.\d+))(?:[EeDd][-+]?\d+)?','match');
values = zeros(1,numel(tokens));
for k = 1:numel(tokens)
    values(k) = str2double(regexprep(tokens{k},'[Dd]','E'));
end
end

function values = fill_constant_missing(values)
finite_index = find(isfinite(values),1,'first');
if isempty(finite_index)
    return;
end
values(~isfinite(values)) = values(finite_index);
end

function values = fill_linear_missing(alpha,values)
valid = isfinite(values);
if sum(valid) >= 2
    values(~valid) = interp1(alpha(valid),values(valid),alpha(~valid), 'linear','extrap');
elseif sum(valid) == 1
    values(~valid) = values(valid);
end
end

function tables = select_best_tables(tables)
if isempty(tables)
    return;
end
keys = strings(numel(tables),1);
for k = 1:numel(tables)
    keys(k) = sprintf('%.5f|%.3f',tables(k).mach,tables(k).altitude_m);
end
unique_keys = unique(keys,'stable');
keep = false(numel(tables),1);
for k = 1:numel(unique_keys)
    candidates = find(keys == unique_keys(k));
    scores = zeros(numel(candidates),1);
    for j = 1:numel(candidates)
        index = candidates(j);
        scores(j) = 1000*tables(index).config_priority+ numel(tables(index).alpha_deg);
    end
    [~,best] = max(scores);
    keep(candidates(best)) = true;
end
tables = tables(keep);
[~,order] = sortrows([[tables.mach].',[tables.altitude_m].'],[1 2]);
tables = tables(order);
end

function reference = validate_reference_geometry(tables)
sref = [tables.sref_m2].';
cbar = [tables.cbar_m].';
bref = [tables.bref_m].';
moment_reference = vertcat(tables.moment_reference_m);
check_consistency(sref,'reference area');
check_consistency(cbar,'reference chord');
check_consistency(bref,'reference span');
for k = 1:3
    check_consistency(moment_reference(:,k), sprintf('moment-reference component %d',k));
end
reference = struct('sref_m2',median(sref), 'cbar_m',median(cbar), 'bref_m',median(bref), 'moment_reference_m',median(moment_reference,1));
end

function check_consistency(values,label)
values = values(isfinite(values));
if isempty(values)
    error('b737datcom:MissingReference','Missing %s.',label);
end
tolerance = max(1e-6,1e-4*max(abs(values)));
if max(values)-min(values) > tolerance
    error('b737datcom:ReferenceMismatch', 'DATCOM files use inconsistent %s values.',label);
end
end

function grid = build_compatibility_grid(static_tables,dynamic_tables)
alpha_deg = unique(vertcat(static_tables.alpha_deg));
mach = unique([static_tables.mach]);
altitude_m = unique([static_tables.altitude_m]);
static_fields = {'CL','CD','CM','CYB','CNB','CLB'};
dynamic_fields = {'CLQ','CMQ','CLAD','CMAD','CLP','CYP','CNP','CNR','CLR'};
size_grid = [numel(alpha_deg),numel(mach),numel(altitude_m)];

grid = struct();
grid.alpha_deg = alpha_deg;
grid.alpha_rad = deg2rad(alpha_deg);
grid.mach_vec = mach(:).';
grid.altitude_m = altitude_m(:).';
grid.alpha_min_deg = NaN(numel(mach),numel(altitude_m));
grid.alpha_max_deg = NaN(numel(mach),numel(altitude_m));
for k = 1:numel(static_fields)
    grid.(static_fields{k}) = NaN(size_grid);
end
for k = 1:numel(dynamic_fields)
    grid.(dynamic_fields{k}) = NaN(size_grid);
end

for k = 1:numel(static_tables)
    table = static_tables(k);
    im = find(abs(mach-table.mach) < 1e-8,1);
    ih = find(abs(altitude_m-table.altitude_m) < 1e-5,1);
    grid.alpha_min_deg(im,ih) = min(table.alpha_deg);
    grid.alpha_max_deg(im,ih) = max(table.alpha_deg);
    for ia = 1:numel(table.alpha_deg)
        ig = find(abs(alpha_deg-table.alpha_deg(ia)) < 1e-8,1);
        for jf = 1:numel(static_fields)
            field = static_fields{jf};
            values = grid.(field);
            source_values = table.(field);
            values(ig,im,ih) = source_values(ia);
            grid.(field) = values;
        end
    end
end

for k = 1:numel(dynamic_tables)
    table = dynamic_tables(k);
    im = find(abs(mach-table.mach) < 1e-8,1);
    ih = find(abs(altitude_m-table.altitude_m) < 1e-5,1);
    if isempty(im) || isempty(ih)
        continue;
    end
    for ia = 1:numel(table.alpha_deg)
        ig = find(abs(alpha_deg-table.alpha_deg(ia)) < 1e-8,1);
        if isempty(ig)
            continue;
        end
        for jf = 1:numel(dynamic_fields)
            field = dynamic_fields{jf};
            values = grid.(field);
            source_values = table.(field);
            values(ig,im,ih) = source_values(ia);
            grid.(field) = values;
        end
    end
end
end

function interpolants = build_interpolants(static_tables,dynamic_tables,altitude_scale)
static_fields = {'CL','CD','CM','CYB','CNB','CLB'};
dynamic_fields = {'CLQ','CMQ','CLAD','CMAD','CLP','CYP','CNP','CNR','CLR'};
interpolants = struct();
for k = 1:numel(static_fields)
    field = static_fields{k};
    interpolants.(field) = make_interpolant(static_tables,field,altitude_scale);
end
for k = 1:numel(dynamic_fields)
    field = dynamic_fields{k};
    interpolants.(field) = make_interpolant(dynamic_tables,field,altitude_scale);
end
end

function pair = make_interpolant(tables,field,altitude_scale)
points = zeros(0,3);
values = zeros(0,1);
for k = 1:numel(tables)
    table = tables(k);
    field_values = table.(field)(:);
    valid = isfinite(table.alpha_deg(:)) & isfinite(field_values);
    count = sum(valid);
    if count == 0
        continue;
    end
    points = [points;[deg2rad(table.alpha_deg(valid)), repmat(table.mach,count,1), repmat(table.altitude_m/altitude_scale,count,1)]]; %#ok<AGROW>
    values = [values;field_values(valid)]; %#ok<AGROW>
end
if size(points,1) < 4
    pair = struct('linear',[],'nearest',[]);
    return;
end
rounded_points = [round(points(:,1),10),round(points(:,2),8),round(points(:,3),10)];
[~,unique_index] = unique(rounded_points,'rows','stable');
points = points(unique_index,:);
values = values(unique_index);
pair = struct();
pair.linear = scatteredInterpolant(points(:,1),points(:,2),points(:,3),values,'linear','none');
pair.nearest = scatteredInterpolant(points(:,1),points(:,2),points(:,3),values,'nearest','nearest');
end

function C = evaluate_lookup(x,u,geom,data)
x = x(:);
u = u(:);
if numel(x) ~= 12 || any(~isfinite(x))
    error('b737datcom:InvalidState','State must be a finite 12-vector.');
end

velocity = x(4:6);
rates = x(10:12);
airspeed = norm(velocity);
alpha = atan2(velocity(3),velocity(1));
beta = atan2(velocity(2),hypot(velocity(1),velocity(3)));
altitude_m = -x(3);
[~,speed_of_sound,~,~] = FlightEnvironment.standard_atmosphere(altitude_m);
Mach = airspeed/max(speed_of_sound,1e-12);
query = [alpha,Mach,altitude_m/data.altitude_scale_m];

[CL,valid_CL] = evaluate_field(data.interpolants.CL,query);
[CD,valid_CD] = evaluate_field(data.interpolants.CD,query);
[Cm,valid_Cm] = evaluate_field(data.interpolants.CM,query);
[CYB,valid_CYB] = evaluate_field(data.interpolants.CYB,query);
[CNB,valid_CNB] = evaluate_field(data.interpolants.CNB,query);
[CLB,valid_CLB] = evaluate_field(data.interpolants.CLB,query);

da = 0;
de = 0;
dr = 0;
if numel(u) >= 1, da = u(1); end
if numel(u) >= 2, de = u(2); end
if numel(u) >= 3, dr = u(3); end

CL = CL+data.control.CLDE*de;
Cm = Cm+data.control.CMDE*de;
CY = CYB*beta+data.control.CYDR*dr;
Cl = CLB*beta+data.control.CLDA*da;
Cn = CNB*beta+data.control.CNDR*dr+data.control.CNDA*da;

dynamic_valid = true;
if airspeed > 1 && any(abs(rates) > 1e-12)
    [~,b,cbar] = geom.get_reference_geometry();
    p_hat = rates(1)*b/(2*airspeed);
    q_hat = rates(2)*cbar/(2*airspeed);
    r_hat = rates(3)*b/(2*airspeed);
    [CLQ,v1] = evaluate_field(data.interpolants.CLQ,query);
    [CMQ,v2] = evaluate_field(data.interpolants.CMQ,query);
    [CLP,v3] = evaluate_field(data.interpolants.CLP,query);
    [CYP,v4] = evaluate_field(data.interpolants.CYP,query);
    [CNP,v5] = evaluate_field(data.interpolants.CNP,query);
    [CNR,v6] = evaluate_field(data.interpolants.CNR,query);
    [CLR,v7] = evaluate_field(data.interpolants.CLR,query);
    CL = CL+CLQ*q_hat;
    Cm = Cm+CMQ*q_hat;
    Cl = Cl+CLP*p_hat+CLR*r_hat;
    CY = CY+CYP*p_hat;
    Cn = Cn+CNP*p_hat+CNR*r_hat;
    dynamic_valid = all([v1,v2,v3,v4,v5,v6,v7]);
end

static_valid = all([valid_CL,valid_CD,valid_Cm, valid_CYB,valid_CNB,valid_CLB]);
C = struct( ...
    'CL',CL,'CD',CD,'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn, ...
    'datcom_valid',static_valid && dynamic_valid, ...
    'datcom_static_valid',static_valid, ...
    'datcom_dynamic_valid',dynamic_valid, ...
    'datcom_extrapolated',~(static_valid && dynamic_valid), ...
    'datcom_query',struct('alpha_rad',alpha,'Mach',Mach, ...
        'altitude_m',altitude_m), ...
    'takeoff_params',data.takeoff_params, ...
    'landing_params',data.landing_params);
end

function [value,valid] = evaluate_field(pair,query)
if isempty(pair.linear)
    value = 0;
    valid = false;
    return;
end
value = pair.linear(query(1),query(2),query(3));
valid = isfinite(value);
if ~valid
    value = pair.nearest(query(1),query(2),query(3));
end
if ~isfinite(value)
    value = 0;
end
end

function control = default_control_model()
control = struct('CLDE',0.00,'CMDE',-1.10, 'CYDR',0.18,'CNDR',-0.070, 'CLDA',0.070,'CNDA',0.010);
end

function control = complete_control_model(control)
defaults = default_control_model();
fields = fieldnames(defaults);
for k = 1:numel(fields)
    field = fields{k};
    if ~isfield(control,field) || isempty(control.(field))
        control.(field) = defaults.(field);
    end
    value = control.(field);
    if ~isscalar(value) || ~isreal(value) || ~isfinite(value)
        error('b737datcom:InvalidControlModel', 'Control derivative %s must be a finite scalar.',field);
    end
end
end

function print_summary(data)
fprintf('\n=== B737 DATCOM DATABASE ===\n');
fprintf('Files               : %d\n',numel(data.source_files));
fprintf('Static tables       : %d\n',numel(data.static_tables));
fprintf('Dynamic tables      : %d\n',numel(data.dynamic_tables));
fprintf('Mach range          : %.3f to %.3f (%d values)\n', min(data.grid.mach_vec),max(data.grid.mach_vec),numel(data.grid.mach_vec));
fprintf('Altitude range      : %.0f to %.0f m (%d values)\n', min(data.grid.altitude_m),max(data.grid.altitude_m), numel(data.grid.altitude_m));
fprintf('Alpha union         : %.1f to %.1f deg (%d values)\n', min(data.grid.alpha_deg),max(data.grid.alpha_deg), numel(data.grid.alpha_deg));
fprintf('Reference S/b/c     : %.3f / %.3f / %.3f m\n', data.reference.sref_m2,data.reference.bref_m,data.reference.cbar_m);
fprintf('Moment reference    : [%.3f %.3f %.3f] m\n', data.reference.moment_reference_m);
fprintf('Control derivatives : %s\n',data.control_model_status);
fprintf('Out-of-range action : %s\n',data.out_of_range_action);
end

function tables = empty_static_tables()
tables = struct('source_file',{},'configuration',{},'case_label',{}, ...
    'config_priority',{},'mach',{},'altitude_m',{},'alpha_deg',{}, ...
    'CD',{},'CL',{},'CM',{},'CYB',{},'CNB',{},'CLB',{}, ...
    'sref_m2',{},'cbar_m',{},'bref_m',{},'moment_reference_m',{});
end

function tables = empty_dynamic_tables()
tables = struct('source_file',{},'configuration',{},'case_label',{}, 'config_priority',{},'mach',{},'altitude_m',{},'alpha_deg',{}, 'CLQ',{},'CMQ',{},'CLAD',{},'CMAD',{},'CLP',{},'CYP',{}, 'CNP',{},'CNR',{},'CLR',{});
end
