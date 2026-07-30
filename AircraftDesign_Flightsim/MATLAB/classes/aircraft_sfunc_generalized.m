function aircraft_sfunc_generalized(block)
% AIRCRAFT_SFUNC_GENERALIZED
% Level-2 MATLAB S-function bridge between Simulink and the Aircraft object.
%
% Inputs:
%   1) x : 12x1 state vector [pos_NED; vel_body; euler; rates_body]
%   2) u : n_total x 1 control vector [control surfaces; propulsive elements]
%
% Outputs:
%   1) F_body : 3x1 total force in body axes [N]
%   2) M_body : 3x1 total moment about CG in body axes [N-m]
%   3) fuel_flow_negative : scalar fuel mass rate [kg/s]
%   4) mass : scalar aircraft mass [kg]
%   5) I_dot : 3x3 inertia derivative [kg-m^2/s]
%   6) I_body : 3x3 inertia matrix about CG/body axes [kg-m^2]
%   7) V_body : 3x1 body velocity [m/s]

setup(block);
end

function setup(block)

block.NumDialogPrms = 0;

n_controls = 4;

if evalin("base", "exist('n_total','var')")
    n_from_base = evalin("base", "n_total");

    if isnumeric(n_from_base) && isscalar(n_from_base) && isfinite(n_from_base) && n_from_base >= 1
        n_controls = double(n_from_base);
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
block.SimStateCompliance = "DefaultSimState";

block.RegBlockMethod("Start",   @Start);
block.RegBlockMethod("Outputs", @Outputs);

end

function Start(block) %#ok<INUSD>

if ~evalin("base", "exist('ac','var')")
    error("aircraft_sfunc_generalized:NoAircraft", "Base workspace variable ac not found. Run setup_f16_simulink_clean first.");
end

ac = evalin("base", "ac");

n_cs    = numel(ac.control_surfaces);
n_pe    = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

assignin("base", "n_cs", n_cs);
assignin("base", "n_pe", n_pe);
assignin("base", "n_total", n_total);

assignin("base", "ac_sfunc", ac);
assignin("base", "n_cs_sfunc", n_cs);
assignin("base", "n_pe_sfunc", n_pe);
assignin("base", "n_total_sfunc", n_total);

end

function Outputs(block)

persistent last_t last_u

if evalin("base", "exist('ac_sfunc','var')")
    ac = evalin("base", "ac_sfunc");
else
    ac = evalin("base", "ac");
end

n_cs    = numel(ac.control_surfaces);
n_pe    = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

t = block.CurrentTime;

% Reset persistent memory at the beginning of a new simulation run.
if isempty(last_t) || t <= 0 || t < last_t
    last_t = t;
    last_u = [];
end

dt = t - last_t;
last_t = t;

%% Read state

x = block.InputPort(1).Data(:);

if numel(x) < 12
    x(end+1:12,1) = 0;
else
    x = x(1:12);
end

x(7) = local_wrap_to_pi(x(7));
x(9) = local_wrap_to_pi(x(9));
x(8) = min(max(x(8), -pi/2 + 0.01), pi/2 - 0.01);

%% Read controls

u_in = block.InputPort(2).Data(:);

if numel(u_in) < n_total
    u_in(end+1:n_total,1) = 0;
else
    u_in = u_in(1:n_total);
end

if isempty(last_u) || numel(last_u) ~= n_total
    last_u = u_in;
end

u_out = u_in;

%% Saturate controls

for i = 1:n_cs
    u_out(i) = min(max(u_out(i), ac.control_surfaces(i).min_deflection), ac.control_surfaces(i).max_deflection);
end

if n_pe > 0
    u_out(n_cs+1:end) = min(max(u_out(n_cs+1:end), 0), 1);
end

%% Optional rate limiting

use_rate_limiting = true;

if evalin("base", "exist('disable_rate_limiting','var')")
    use_rate_limiting = ~logical(evalin("base", "disable_rate_limiting"));
end

if use_rate_limiting
    max_du = zeros(size(u_out));

    if n_cs > 0
        max_du(1:n_cs) = deg2rad(60) * dt;
    end

    if n_pe > 0
        max_du(n_cs+1:end) = 0.5 * dt;
    end

    du = min(max(u_out - last_u, -max_du), max_du);
    u_out = last_u + du;
end

last_u = u_out;

%% Update aircraft

ac.state.set_full_state(x);
ac.set_controls_from_vector(u_out);

%% Mass and CG

[m_total, cg_body, I_body] = ac.compute_total_mass_properties(x);

if ac.has_frame("cg")
    ac.update_frame_position("cg", cg_body);
else
    ac.add_frame("cg", ac.body_frame_name, cg_body, @(x) eye(3));
end

if ac.has_frame("gravity_cg")
    ac.update_frame_position("gravity_cg", cg_body);
end

%% Total loads about CG

old_ref = ac.reference_frame_name;

try
    ac.set_reference_frame("cg");

    [F_body, M_body, fuel_flow] = ac.compute_total_loads(x, u_out);

    ac.set_reference_frame(old_ref);

catch ME
    try
        ac.set_reference_frame(old_ref);
    catch
    end

    error("aircraft_sfunc_generalized:ForceMomentFailure", "Failed to compute total loads at t = %.6f s.\n%s", t, ME.message);
end

if isempty(F_body) || numel(F_body) ~= 3
    error("aircraft_sfunc_generalized:BadForceSize", "F_body must be 3x1.");
end

if isempty(M_body) || numel(M_body) ~= 3
    error("aircraft_sfunc_generalized:BadMomentSize", "M_body must be 3x1.");
end

if isempty(fuel_flow) || ~isfinite(fuel_flow)
    fuel_flow = 0;
end

I_dot = zeros(3,3);

%% Outputs

block.OutputPort(1).Data = F_body(:);
block.OutputPort(2).Data = M_body(:);
block.OutputPort(3).Data = -fuel_flow;
block.OutputPort(4).Data = m_total;
block.OutputPort(5).Data = I_dot;
block.OutputPort(6).Data = I_body;
block.OutputPort(7).Data = x(4:6);

%% Debug exports

assignin("base", "sfunc_last_x", x);
assignin("base", "sfunc_last_u", u_out);
assignin("base", "sfunc_last_F_body", F_body(:));
assignin("base", "sfunc_last_M_body", M_body(:));
assignin("base", "sfunc_last_mass", m_total);
assignin("base", "sfunc_last_I_body", I_body);
assignin("base", "sfunc_last_cg_body", cg_body);

%% Ground warning

altitude = -x(3);

if altitude < -5
    warning("aircraft_sfunc_generalized:GroundCollision", "Ground collision warning at t = %.2f s, altitude = %.1f m", t, altitude);
end

end

function y = local_wrap_to_pi(x)

y = mod(x + pi, 2*pi) - pi;

if y == -pi
    y = pi;
end

end