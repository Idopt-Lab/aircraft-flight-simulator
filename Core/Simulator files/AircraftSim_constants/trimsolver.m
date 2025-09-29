clear; clc;

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

function [trim_state, trim_controls, converged] = find_trim_condition(aircraft, altitude, mach_number)
    [T, a, P, rho] = atmosisa(altitude);
    V = mach_number * a;
    
    options = optimset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIter', 1000);
    
    x0 = [deg2rad(3.0); deg2rad(-1.0); 0.5];
    
    [x_trim, fval, exitflag] = fsolve(@(x) trim_equations(x, aircraft, altitude, V), x0, options);
    
    converged = (exitflag > 0) && (norm(fval) < 1e-3);
    
    if converged
        alpha_trim = x_trim(1);
        elevator_trim = x_trim(2);
        throttle_trim = x_trim(3);
        
        trim_state = zeros(12,1);
        trim_state(1) = 0; trim_state(2) = 0; trim_state(3) = -altitude;
        trim_state(4) = V*cos(alpha_trim); trim_state(5) = 0; trim_state(6) = V*sin(alpha_trim);
        trim_state(7) = 0; trim_state(8) = alpha_trim; trim_state(9) = 0;
        trim_state(10) = 0; trim_state(11) = 0; trim_state(12) = 0;
        
        trim_controls = struct();
        trim_controls.alpha_deg = rad2deg(alpha_trim);
        trim_controls.elevator_deg = rad2deg(elevator_trim);
        trim_controls.throttle = throttle_trim;
        trim_controls.aileron_deg = 0;
        trim_controls.rudder_deg = 0;
    else
        trim_state = [];
        trim_controls = [];
    end
end

function residuals = trim_equations(x, aircraft, altitude, V)
    alpha = x(1);
    elevator_deflection = x(2);
    throttle = x(3);
    
    throttle = max(0, min(1, throttle));
    alpha = max(deg2rad(-5), min(deg2rad(20), alpha));
    elevator_deflection = max(deg2rad(-25), min(deg2rad(25), elevator_deflection));
    
    state = zeros(12,1);
    state(1) = 0; state(2) = 0; state(3) = -altitude;
    state(4) = V*cos(alpha); state(5) = 0; state(6) = V*sin(alpha);
    state(7) = 0; state(8) = alpha; state(9) = 0;
    state(10) = 0; state(11) = 0; state(12) = 0;
    aircraft.state.set_full_state(state);
    
    aircraft.set_control_by_name('aileron', 0);
    aircraft.set_control_by_name('elevator', elevator_deflection);
    aircraft.set_control_by_name('rudder', 0);
    aircraft.set_control_by_name('F100_PW_200', throttle);
    
    state_current = aircraft.state.get_full_state();
    controls_current = aircraft.control.get_full_controls();
    [F_aero, M_aero] = aircraft.aero.calculate_forces_moments(state_current, controls_current, aircraft.geometry, aircraft.control_surfaces);
    [F_prop, M_prop] = PropulsiveElement.calculate_total_forces(aircraft.propulsive_elements, state_current, controls_current);
    [F_weight, M_weight] = aircraft.mass.calculate_weight_forces(state_current, 'constant');
    
    F_total = F_aero + F_prop + F_weight;
    M_total = M_aero + M_prop + M_weight;
    
    residuals(1) = F_total(1);
    residuals(2) = F_total(3);
    residuals(3) = M_total(2);
end

alt_m = 9144;
M = 0.60;

[trim_state, trim_controls, converged] = find_trim_condition(f16, alt_m, M);

if converged
    f16.state.set_full_state(trim_state);
    f16.set_control_by_name('aileron', deg2rad(trim_controls.aileron_deg));
    f16.set_control_by_name('elevator', deg2rad(trim_controls.elevator_deg));
    f16.set_control_by_name('rudder', deg2rad(trim_controls.rudder_deg));
    f16.set_control_by_name('F100_PW_200', trim_controls.throttle);
    
    state = f16.state.get_full_state();
    controls = f16.control.get_full_controls();
    [F_aero, M_aero] = f16.aero.calculate_forces_moments(state, controls, f16.geometry, f16.control_surfaces);
    [F_prop, M_prop] = PropulsiveElement.calculate_total_forces(f16.propulsive_elements, state, controls);
    [F_weight, M_weight] = f16.mass.calculate_weight_forces(state, 'constant');
    
    F_total = F_aero + F_prop + F_weight;
    M_total = M_aero + M_prop + M_weight;
    
    Initialpos = [f16.state.state(1), f16.state.state(2), f16.state.state(3)];
    InitialVel = [f16.state.state(4), f16.state.state(5), f16.state.state(6)];
    InitialOri = [f16.state.state(7), f16.state.state(8), f16.state.state(9)];
    InitialRot = [f16.state.state(10), f16.state.state(11), f16.state.state(12)];
    
    trim_state
    trim_controls
    F_total
    M_total
else
    fprintf('Trim solution did not converge\n');
end