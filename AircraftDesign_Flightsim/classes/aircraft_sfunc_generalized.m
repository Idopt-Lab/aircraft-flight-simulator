function aircraft_sfunc_generalized(block)
setup(block);
end

function setup(block)
block.NumDialogPrms = 0;

n_controls = 4;
if evalin('base','exist(''n_total'',''var'')')
    v = evalin('base','n_total');
    if isnumeric(v) && isscalar(v) && isfinite(v) && v >= 1
        n_controls = double(v);
    end
end

block.NumInputPorts  = 2;
block.NumOutputPorts = 7;

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions        = 12;
block.InputPort(1).DirectFeedthrough = true;

block.InputPort(2).Dimensions        = n_controls;
block.InputPort(2).DirectFeedthrough = true;

block.OutputPort(1).Dimensions = 3;
block.OutputPort(2).Dimensions = 3;
block.OutputPort(3).Dimensions = 1;
block.OutputPort(4).Dimensions = 1;
block.OutputPort(5).Dimensions = [3 3];
block.OutputPort(6).Dimensions = [3 3];
block.OutputPort(7).Dimensions = 3;

block.SampleTimes = [-1 0];
block.SimStateCompliance = 'DefaultSimState';

block.RegBlockMethod('Start',   @Start);
block.RegBlockMethod('Outputs', @Outputs);
end

function Start(block)
persistent ac autopilot n_total n_cs n_pe

if isempty(ac)
    if evalin('base','exist(''ac'',''var'')')
        ac = evalin('base','ac');
    else
        error('aircraft_sfunc_generalized:NoAircraft','ac not found');
    end

    if evalin('base','exist(''autopilot'',''var'')')
        autopilot = evalin('base','autopilot');
    else
        autopilot = [];
    end

    n_cs = numel(ac.control_surfaces);
    n_pe = numel(ac.propulsive_elements);
    n_total = n_cs + n_pe;

    assignin('base','n_total',n_total);
end

assignin('base','ac_sfunc',ac);
assignin('base','autopilot_sfunc',autopilot);
assignin('base','n_total_sfunc',n_total);
assignin('base','n_cs_sfunc',n_cs);
assignin('base','n_pe_sfunc',n_pe);
end

function Outputs(block)
persistent ac autopilot n_total n_cs n_pe last_t last_u ap_prev_enabled

if isempty(last_t), last_t = block.CurrentTime; end
if isempty(ap_prev_enabled), ap_prev_enabled = false; end

if isempty(ac)
    ac = evalin('base','ac_sfunc');
    autopilot = [];
    if evalin('base','exist(''autopilot_sfunc'',''var'')')
        autopilot = evalin('base','autopilot_sfunc');
    end
    n_total = evalin('base','n_total_sfunc');
    n_cs    = evalin('base','n_cs_sfunc');
    n_pe    = evalin('base','n_pe_sfunc');
end

t  = block.CurrentTime;
dt = t - last_t;
if ~isfinite(dt) || dt <= 0, dt = 0.01; end
dt = min(max(dt,1e-4), 0.05);
last_t = t;

x = block.InputPort(1).Data(:);
if numel(x) < 12
    x(end+1:12,1) = 0;
else
    x = x(1:12);
end

x(7) = wrapToPi(x(7));
x(9) = wrapToPi(x(9));
x(8) = min(max(x(8), -pi/2+0.01), pi/2-0.01);

u_in = block.InputPort(2).Data(:);
if numel(u_in) < n_total
    u_in(end+1:n_total,1) = 0;
else
    u_in = u_in(1:n_total);
end

if isempty(last_u) || numel(last_u) ~= n_total
    last_u = u_in;
end

altitude = -x(3);

ap_en = 0;
ap_md = "off";
ap_t_on = 0;
ap_min_alt = 0;

if evalin('base','exist(''autopilot_enabled'',''var'')'), ap_en = double(evalin('base','autopilot_enabled')); end
if evalin('base','exist(''autopilot_mode'',''var'')'), ap_md = string(evalin('base','autopilot_mode')); end
if evalin('base','exist(''autopilot_enable_time'',''var'')'), ap_t_on = double(evalin('base','autopilot_enable_time')); end
if evalin('base','exist(''autopilot_min_alt'',''var'')'), ap_min_alt = double(evalin('base','autopilot_min_alt')); end

ap_active = logical(ap_en) && (t >= ap_t_on) && (altitude >= ap_min_alt) && (ap_md ~= "off");

u_cmd = u_in;

if ~isempty(autopilot) && isobject(autopilot) && isprop(autopilot,'enabled') && isprop(autopilot,'mode')
    autopilot.mode = ap_md;
    autopilot.enabled = ap_active;

    if autopilot.enabled && autopilot.mode ~= "off"
        if evalin('base','exist(''mission'',''var'')')
            mission = evalin('base','mission');
            t_m = min(max(t, 0), mission.total_duration);
            autopilot.target_altitude = interp1(mission.time_vector, mission.altitude_profile, t_m, 'linear', 'extrap');
            autopilot.target_speed = interp1(mission.time_vector, mission.velocity_profile, t_m, 'linear', 'extrap');
        end

        if ~ap_prev_enabled
            try
                autopilot.initialize_from_state(x, u_in);
            catch
            end
        end

        try
            u_cmd = autopilot.compute_control(x, u_in, dt);
        catch
            u_cmd = u_in;
        end
    end
end

ap_prev_enabled = ap_active;

u_out = u_cmd;

for i = 1:n_cs
    mn = ac.control_surfaces(i).min_deflection;
    mx = ac.control_surfaces(i).max_deflection;
    u_out(i) = min(max(u_out(i), mn), mx);
end

if n_pe > 0
    idx = (n_cs+1):(n_cs+n_pe);
    u_out(idx) = min(max(u_out(idx), 0), 1);
end

if ~isempty(last_u)
    max_du = zeros(size(u_out));
    max_du(1:n_cs) = deg2rad(60) * dt;
    if n_pe > 0
        max_du(n_cs+1:end) = 0.5 * dt;
    end
    du = u_out - last_u;
    du = min(max(du, -max_du), max_du);
    u_out = last_u + du;
end

last_u = u_out;

ac.state.set_full_state(x);
ac.set_controls_from_vector(u_out);

[F_total, M_total, fuel_flow] = ac.calculate_total_forces_moments_with_gravity();

phi = x(7); th = x(8); ps = x(9);

cph = cos(phi); sph = sin(phi);
cth = cos(th);  sth = sin(th);
cps = cos(ps);  sps = sin(ps);

Cbn = [ cth*cps,               cth*sps,              -sth;
        sph*sth*cps-cph*sps,   sph*sth*sps+cph*cps,  sph*cth;
        cph*sth*cps+sph*sps,   cph*sth*sps-sph*cps,  cph*cth ];

Vb = x(4:6);
Vned = Cbn.' * Vb;
Vd = Vned(3);

k_g = 0;
c_g = 0;
if evalin('base','exist(''ground_k'',''var'')'), k_g = double(evalin('base','ground_k')); end
if evalin('base','exist(''ground_c'',''var'')'), c_g = double(evalin('base','ground_c')); end

if x(3) > 0
    depth = x(3);
    Fn_ned = [0;0;-(k_g*depth + c_g*max(Vd,0))];
    Fg_b = Cbn * Fn_ned;
    F_total = F_total + Fg_b;
end

m     = ac.mass.get_total_mass();
I_mat = ac.mass.get_inertia_matrix();
I_dot = zeros(3,3);

block.OutputPort(1).Data = F_total(:);
block.OutputPort(2).Data = M_total(:);
block.OutputPort(3).Data = -fuel_flow;
block.OutputPort(4).Data = m;
block.OutputPort(5).Data = I_dot;
block.OutputPort(6).Data = I_mat;
block.OutputPort(7).Data = Vb(:);

altitude = -x(3);
if altitude < -5
    warning('aircraft_sfunc:GroundCollision','GROUND COLLISION at t=%.2f s, altitude=%.1f m', t, altitude);
end
end
