classdef (Abstract) Mass < handle
    % MASS  Abstract base class for aircraft mass modeling.
    %
    %   Defines the common interface for aircraft mass models. Each subclass
    %   must provide total mass, empty mass, current center of gravity, and
    %   inertia tensor information. This base class also provides shared fuel
    %   management, gravity force/moment computation, and a protected parallel
    %   axis theorem utility.
    %
    %   Coordinate and modeling assumptions:
    %     1. Body axes follow the standard aircraft convention:
    %        x forward, y right, z downward.
    %     2. The inertial/navigation frame is assumed to be NED.
    %     3. All CG locations, fuel tank locations, and reference points are
    %        defined in body axes [m].
    %     4. Euler angles are ordered as [phi; theta; psi] = roll, pitch, yaw [rad].
    %     5. Gravity is modeled as a constant acceleration acting in the NED
    %        down direction.
    %     6. Gravity acts through the aircraft CG. If moments are evaluated about
    %        a reference point other than the CG, the gravity moment is computed as
    %        M_g = r_cg x F_g.
    %     7. Fuel is treated as lumped mass at each tank location unless a
    %        subclass provides a more detailed fuel inertia model.
    %
    %   Subclasses:
    %     SimpleMass     - Top-down mass model using specified mass, CG, and inertia.
    %     ComponentMass  - Bottom-up mass model built from discrete components.
    %
    %   References:
    %     Stevens, B. L., Lewis, F. L., and Johnson, E. N.,
    %     Aircraft Control and Simulation, 3rd ed., Wiley, 2015.
    %
    %     Etkin, B. and Reid, L. D.,
    %     Dynamics of Flight: Stability and Control, 3rd ed., Wiley, 1996.
    %
    %   See also: SimpleMass, ComponentMass, Aircraft
    properties
        % Current total fuel mass [kg]
        fuel_mass = 0

        % Total fuel capacity across all tanks [kg]
        max_fuel = 0

        % Maximum take-off weight [kg] (informational, not enforced)
        mtow = 0

        % Nominal fuel burn rate [kg/s] (informational)
        fuel_burn_rate = 0

        % Struct array of fuel tanks: location, capacity, current, type
        fuel_tanks = []

        % Current CG location in body axes [m] (auto-updated)
        cg_current = [0 0 0]

        % Current inertia tensor in body axes [kg-m^2] (auto-updated)
        I_body_current = zeros(3,3)

        % Gravitational acceleration [m/s^2]
        g = 9.80665
    end

    methods (Abstract)
        % Get total aircraft mass (empty + fuel) [kg]
        m = get_total_mass(obj)

        % Get empty aircraft mass [kg]
        m = get_empty_mass(obj)

        % Get current CG location in body axes [m]
        cg = get_cg(obj)

        % Get current inertia matrix [kg-m^2]
        I = get_inertia_matrix(obj)
    end

    methods (Abstract, Access = protected)
        % Recompute CG and inertia for current configuration
        update_current_properties(obj)
    end

    methods
        % ── Fuel management (common to all subclasses) ───────────────────

        function add_fuel_tank(obj, location, capacity, tank_type)
            % ADD_FUEL_TANK  Register a fuel tank.
            %
            %   Inputs:
            %     location  - 1x3 body-frame position [m]
            %     capacity  - maximum fuel mass [kg]
            %     tank_type - 'distributed' (default) or 'sphere'

            if nargin < 4, tank_type = 'distributed'; end
            n = numel(obj.fuel_tanks);
            obj.fuel_tanks(n+1).location = location(:)';
            obj.fuel_tanks(n+1).capacity = capacity;
            obj.fuel_tanks(n+1).current  = 0;
            obj.fuel_tanks(n+1).type     = tank_type;
            obj.max_fuel = obj.max_fuel + capacity;
            obj.update_current_properties();
        end

        function set_fuel_properties(obj, max_fuel, mtow)
            % SET_FUEL_PROPERTIES  Set total fuel capacity and MTOW.

            obj.max_fuel = max_fuel;
            obj.mtow     = mtow;
            obj.update_current_properties();
        end

        function set_fuel_capacity(obj, max_fuel, mtow)
            % SET_FUEL_CAPACITY  Alias for set_fuel_properties.

            obj.max_fuel = max_fuel;
            if nargin > 2 && ~isempty(mtow)
                obj.mtow = mtow;
            end
            obj.update_current_properties();
        end

        function set_fuel_mass(obj, m_fuel)
            % SET_FUEL_MASS  Set total fuel and distribute across tanks sequentially.

            obj.fuel_mass = max(0, min(obj.max_fuel, m_fuel));
            if ~isempty(obj.fuel_tanks)
                remaining = obj.fuel_mass;
                for i = 1:numel(obj.fuel_tanks)
                    fill = min(remaining, obj.fuel_tanks(i).capacity);
                    obj.fuel_tanks(i).current = fill;
                    remaining = remaining - fill;
                end
            end
            obj.update_current_properties();
        end

        function set_fuel_burn_rate(obj, rate)
            % SET_FUEL_BURN_RATE  Store nominal fuel burn rate (informational).

            obj.fuel_burn_rate = rate;
        end

        function update_fuel_mass(obj, dm)
            % UPDATE_FUEL_MASS  Add or remove fuel and update tank states.
            %
            %   Input:
            %     dm - fuel mass change [kg] (negative = burn, positive = fill)

            obj.fuel_mass = max(0, obj.fuel_mass + dm);
            if ~isempty(obj.fuel_tanks)
                if dm < 0
                    burn = -dm;
                    for i = 1:numel(obj.fuel_tanks)
                        if burn <= 0, break; end
                        take = min(burn, obj.fuel_tanks(i).current);
                        obj.fuel_tanks(i).current = obj.fuel_tanks(i).current - take;
                        burn = burn - take;
                    end
                else
                    fill = dm;
                    for i = 1:numel(obj.fuel_tanks)
                        if fill <= 0, break; end
                        space = obj.fuel_tanks(i).capacity - obj.fuel_tanks(i).current;
                        add   = min(fill, space);
                        obj.fuel_tanks(i).current = obj.fuel_tanks(i).current + add;
                        fill = fill - add;
                    end
                end
            end
            obj.update_current_properties();
        end

        function update_fuel_from_flow(obj, fuel_flow, dt)
            % UPDATE_FUEL_FROM_FLOW  Burn fuel based on flow rate and time step.

            obj.update_fuel_mass(-fuel_flow * dt);
        end

        function m = get_fuel_mass(obj)
            % GET_FUEL_MASS  Return current fuel mass [kg].

            m = obj.fuel_mass;
        end

    

        function [F_g, M_g] = get_gravity_force_moment(obj, euler_angles, ref_point)
            % GET_GRAVITY_FORCE_MOMENT  Compute gravity force and moment in body axes.
            %
            %   Gravity force is projected from NED to body axes using Euler angles.
            %   Gravity moment is computed about a specified reference point using
            %   M_g = r_cg × F_g, where r_cg is the vector from ref_point to CG.
            %
            %   Inputs:
            %     euler_angles - 3x1 [phi; theta; psi] [rad]
            %     ref_point    - 3x1 or 1x3 reference point in body axes [m]
            %                    (optional, default [0; 0; 0] = nose/origin)
            %
            %   Outputs:
            %     F_g - 3x1 gravity force in body axes [N]
            %     M_g - 3x1 gravity moment about ref_point in body axes [N-m]
            %
            %   %   References:
            %     The body-axis gravity projection follows the standard 3-2-1 Euler-angle
            %     aircraft equations of motion for a flat-Earth NED frame.
            %     See Stevens, Lewis, and Johnson, Aircraft Control and Simulation, 3rd ed.

            if nargin < 3 || isempty(ref_point)
                ref_point = [0; 0; 0];
            end

            euler_angles = euler_angles(:);
            ref_point    = ref_point(:);

            m     = obj.get_total_mass();
            phi   = euler_angles(1);
            theta = euler_angles(2);

            % Gravity force in body axes (psi cancels for flat Earth)
            F_g = m * obj.g * [-sin(theta);
                sin(phi) * cos(theta);
                cos(phi) * cos(theta)];

            % Moment arm from reference point to CG
            r_cg = obj.cg_current(:) - ref_point;

            % Gravity moment about reference point: M_g = r_cg × F_g
            M_g = cross(r_cg, F_g);
        end

        function set_gravity(obj, g_value)
            % SET_GRAVITY  Override gravitational acceleration.
            %
            %   Input:
            %     g_value - gravitational acceleration [m/s^2]
            %               (default 9.80665 m/s^2 at sea level)

            obj.g = g_value;
        end

    end

    methods (Access = protected)

        function I_shifted = parallel_axis_theorem(~, I_cm, mass, dcg)
            % PARALLEL_AXIS_THEOREM  Shift inertia tensor from component CG to new reference.
            %
            %   I_new = I_cm + m*([d^2*I - d*d^T])
            %
            %   Inputs:
            %     I_cm - 3x3 inertia at component CG [kg-m^2]
            %     mass - component mass [kg]
            %     dcg  - 3x1 displacement vector [m]
            %
            %   Output:
            %     I_shifted - 3x3 inertia at new reference [kg-m^2]

            dcg = dcg(:);
            dx = dcg(1); dy = dcg(2); dz = dcg(3);
            d2 = dx^2 + dy^2 + dz^2;

            I_parallel = mass * [d2-dx^2, -dx*dy,  -dx*dz;
                -dx*dy,  d2-dy^2, -dy*dz;
                -dx*dz,  -dy*dz,  d2-dz^2];

            I_shifted = I_cm + I_parallel;
        end

    end
end
