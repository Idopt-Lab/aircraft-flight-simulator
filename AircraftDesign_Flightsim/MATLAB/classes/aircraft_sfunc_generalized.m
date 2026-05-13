function aircraft_sfunc_generalized(block)
% AIRCRAFT_SFUNC_GENERALIZED  Simulink Level-2 S-function for the aircraft model.
%
%   Bridges the Aircraft object (in the MATLAB base workspace) to a
%   Simulink simulation. At each time step it reads the state and control
%   inputs, evaluates total forces and moments via the Aircraft class, and
%   outputs the results for the equations-of-motion integrator.
%
%   Base workspace variables required before simulation:
%     ac       - Aircraft object (fully configured)
%     n_total  - total number of controls (n_cs + n_pe)
%
%   
%
%   See also: Aircraft, setup, Start, Outputs

setup(block);
end

function setup(block)
% SETUP  Register block ports, dimensions, and sample time.
%
%   Input ports:
%     1 - state vector x, 12x1 [pos(3); vel(3); euler(3); rates(3)]
%     2 - control vector u, n_total x 1
%
%   Output ports:
%     1 - total force vector F_total, 3x1 [N]
%     2 - total moment vector M_ext, 3x1 [N-m]
%     3 - fuel flow rate (negative), scalar [kg/s]
%     4 - total aircraft mass, scalar [kg]
%     5 - inertia rate matrix I_dot, 3x3 [kg-m^2/s]
%     6 - inertia matrix I, 3x3 [kg-m^2]
%     7 - body velocity vector, 3x1 [m/s]

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
% START  Initialise persistent aircraft object at simulation start.
%
%   Loads the Aircraft object from the base workspace, extracts control
%   surface and propulsion element counts, and re-exports them to the
%   base workspace for use by other blocks.

persistent ac n_total n_cs n_pe

if isempty(ac)
    if evalin('base','exist(''ac'',''var'')')
        ac = evalin('base','ac');
    else
        error('aircraft_sfunc_generalized:NoAircraft','ac not found');
    end

    n_cs    = numel(ac.control_surfaces);
    n_pe    = numel(ac.propulsive_elements);
    n_total = n_cs + n_pe;

    assignin('base','n_total', n_total);
    assignin('base','n_cs',    n_cs);
    assignin('base','n_pe',    n_pe);
end

assignin('base','ac_sfunc',      ac);
assignin('base','n_total_sfunc', n_total);
assignin('base','n_cs_sfunc',    n_cs);
assignin('base','n_pe_sfunc',    n_pe);
end

function Outputs(block)
% OUTPUTS  Compute and write output ports at each simulation time step.
%
%   At each step this function:
%     1. Reads state x and control u from input ports.
%     2. Wraps Euler angles to valid ranges.
%     3. Saturates controls to surface deflection and throttle limits.
%     4. Applies rate limiting to control inputs (60 deg/s for surfaces,
%        0.5/s for throttle) unless disable_rate_limiting is set.
%     5. Updates the Aircraft object state and controls.
%     6. Evaluates total forces and moments including ground contact.
%     7. Writes forces, moments, fuel flow, mass, inertia, and velocity
%        to the seven output ports.
%     8. Issues a warning if the aircraft altitude goes below -5 m.

persistent ac n_total n_cs n_pe last_t last_u

if isempty(last_t), last_t = block.CurrentTime; end

if isempty(ac)
    ac      = evalin('base','ac_sfunc');
    n_total = evalin('base','n_total_sfunc');
    n_cs    = evalin('base','n_cs_sfunc');
    n_pe    = evalin('base','n_pe_sfunc');
end

t  = block.CurrentTime;
dt = t - last_t;
if ~isfinite(dt) || dt <= 0, dt = 0.01; end
dt = min(max(dt, 1e-4), 0.05);
last_t = t;

x = block.InputPort(1).Data(:);
if numel(x) < 12, x(end+1:12,1) = 0; else, x = x(1:12); end

x(7) = wrapToPi(x(7));
x(9) = wrapToPi(x(9));
x(8) = min(max(x(8), -pi/2+0.01), pi/2-0.01);

u_in = block.InputPort(2).Data(:);
if numel(u_in) < n_total, u_in(end+1:n_total,1) = 0; else, u_in = u_in(1:n_total); end

if isempty(last_u) || numel(last_u) ~= n_total
    last_u = u_in;
end

u_out = u_in;

for i = 1:n_cs
    u_out(i) = min(max(u_out(i), ac.control_surfaces(i).min_deflection), ...
                                  ac.control_surfaces(i).max_deflection);
end
if n_pe > 0
    u_out(n_cs+1:end) = min(max(u_out(n_cs+1:end), 0), 1);
end

use_rate_limiting = 1;
if evalin('base','exist(''disable_rate_limiting'',''var'')')
    use_rate_limiting = ~double(evalin('base','disable_rate_limiting'));
end

if use_rate_limiting
    max_du = zeros(size(u_out));
    max_du(1:n_cs) = deg2rad(60) * dt;
    if n_pe > 0, max_du(n_cs+1:end) = 0.5 * dt; end
    du    = min(max(u_out - last_u, -max_du), max_du);
    u_out = last_u + du;
end

last_u = u_out;

ac.state.set_full_state(x);
ac.set_controls_from_vector(u_out);

k_g = 0; c_g = 0;
if evalin('base','exist(''ground_k'',''var'')'), k_g = double(evalin('base','ground_k')); end
if evalin('base','exist(''ground_c'',''var'')'), c_g = double(evalin('base','ground_c')); end

[F_total, M_ext, fuel_flow, m] = ac.calculate_total_forces_moments_with_ground(k_g, c_g);

I_mat = ac.mass.get_inertia_matrix();
I_dot = zeros(3,3);

block.OutputPort(1).Data = F_total(:);
block.OutputPort(2).Data = M_ext(:);
block.OutputPort(3).Data = -fuel_flow;
block.OutputPort(4).Data = m;
block.OutputPort(5).Data = I_dot;
block.OutputPort(6).Data = I_mat;
block.OutputPort(7).Data = x(4:6);

if -x(3) < -5
    warning('aircraft_sfunc:GroundCollision', ...
        'GROUND COLLISION at t=%.2f s, altitude=%.1f m', t, -x(3));
end
end