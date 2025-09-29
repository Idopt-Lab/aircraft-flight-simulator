function y = aircraft_interface_sim(u)

global AIRCRAFT_OBJ

state_vec = u(1:12);
control_vec = u(13:16);

if ~isempty(AIRCRAFT_OBJ)
    [F_aero, M_aero] = AIRCRAFT_OBJ.aero.calculate_forces_moments(state_vec, control_vec, AIRCRAFT_OBJ.geometry, AIRCRAFT_OBJ.control_surfaces);
    [F_prop, M_prop] = PropulsiveElement.calculate_total_forces(AIRCRAFT_OBJ.propulsive_elements, state_vec, control_vec);
    [F_weight, M_weight] = AIRCRAFT_OBJ.mass.calculate_weight_forces(state_vec, 'constant');
    
    F_total = F_aero + F_prop + F_weight;
    M_total = M_aero + M_prop + M_weight;
else
    F_total = [0; 0; 0];
    M_total = [0; 0; 0];
end

y = [F_total; M_total];

end