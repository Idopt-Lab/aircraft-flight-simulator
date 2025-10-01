clear; clc;
global AIRCRAFT_OBJ

f16 = Aircraft();

weight_N = 91188;
mass_kg = weight_N / 9.81;
Ixx = 12875; Iyy = 75674; Izz = 85552; Ixz = 1331; Iyz = 0; Ixy = 0;
cg = [0, 0, 0];
f16.mass.set_mass_properties(mass_kg, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, cg);

S = 27.87; b = 9.144; c = 3.45;
f16.geometry.add_component('wing', 'main_wing', S, b, c, true, [0, 0, 0]);

aileron = ControlSurface('aileron');
aileron.set_properties('aileron', [0.318 0 -0.063], deg2rad(20), deg2rad(-20), [1 0 0], 'aileron');
f16.add_control_surface(aileron);

elevator = ControlSurface('elevator');
elevator.set_properties('elevator', [0 -1.315 0], deg2rad(25), deg2rad(-25), [0 1 0], 'elevator');
f16.add_control_surface(elevator);

rudder = ControlSurface('rudder');
rudder.set_properties('rudder', [0.064 0 -0.084], deg2rad(30), deg2rad(-30), [0 0 1], 'rudder');
f16.add_control_surface(rudder);

f100 = PropulsiveElement('engine');
f100.set_properties('F100_PW_200', 23830, [-1.8 0 0], [1 0 0], 0.477);
f16.add_propulsive_element(f100);

lookup = F16Lookup();
f16.aero.set_lookup(lookup);

alt_m = 9144;
M = 0.60;
alpha_deg = 5.2402;
[T, a, P, rho] = atmosisa(alt_m);
V = M * a;
alpha = deg2rad(alpha_deg);

initial_state = zeros(12,1);
initial_state(1) = 0; initial_state(2) = 0; initial_state(3) = -alt_m;
initial_state(4) = V*cos(alpha); initial_state(5) = 0; initial_state(6) = V*sin(alpha);
initial_state(7) = 0; initial_state(8) = alpha; initial_state(9) = 0;
initial_state(10) = 0; initial_state(11) = 0; initial_state(12) = 0;

initial_controls = zeros(4,1);
initial_controls(1) = deg2rad(0.0);
initial_controls(2) = deg2rad(0.6575);
initial_controls(3) = deg2rad(0.0);
initial_controls(4) = 0.3223;

f16.state.set_full_state(initial_state);
f16.control.set_control(1, initial_controls(1));
f16.control.set_control(2, initial_controls(2));
f16.control.set_control(3, initial_controls(3));
f16.control.set_control(4, initial_controls(4));

AIRCRAFT_OBJ = f16;

combined_input = [initial_state; initial_controls];
try
    result = aircraft_interface_sim(combined_input);
    F = result(1:3);
    M = result(4:6);
    fprintf('Interface function test successful.\n');
    fprintf('Forces: [%.3f %.3f %.3f] N\n', F(1), F(2), F(3));
    fprintf('Moments: [%.3f %.3f %.3f] N-m\n', M(1), M(2), M(3));
catch ME
    error('Interface function failed: %s', ME.message);
end

Initialpos = initial_state(1:3);
InitialVel = initial_state(4:6);
InitialOri = initial_state(7:9);
InitialRot = initial_state(10:12);
Inertia = f16.mass.get_inertia_matrix();
total_mass = f16.mass.get_total_mass();

assignin('base', 'Initialpos', Initialpos);
assignin('base', 'InitialVel', InitialVel);
assignin('base', 'InitialOri', InitialOri);
assignin('base', 'InitialRot', InitialRot);
assignin('base', 'Inertia', Inertia);
assignin('base', 'mass', total_mass);

model_name = 'Flightsim2';

try
    if ~bdIsLoaded(model_name)
        load_system(model_name);
    end
    
    sixdof_blocks = [
        find_system(model_name, 'BlockType', 'S-Function');
        find_system(model_name, 'MaskType', '6DOF (Euler Angles)');
        find_system(model_name, 'ReferenceBlock', 'aeroblkspaceenv/6DOF (Euler Angles)')
    ];
    
    for i = 1:length(sixdof_blocks)
        block = sixdof_blocks{i};
        try
            params = get_param(block, 'DialogParameters');
            if isfield(params, 'rp_0')
                set_param(block, 'rp_0', 'Initialpos');
                set_param(block, 'vb_0', 'InitialVel');
                set_param(block, 'eul_0', 'InitialOri');
                set_param(block, 'wb_0', 'InitialRot');
                fprintf('6DOF configured: %s\n', block);
                break;
            end
        catch
            continue;
        end
    end
catch
    fprintf('Auto-config failed\n');
end

set_param(model_name, 'StopTime', '100');
set_param(model_name, 'FixedStep', '0.01');
try
    sim(model_name);
    fprintf('Simulation completed\n');
catch ME
    fprintf('Simulation error: %s\n', ME.message);
end