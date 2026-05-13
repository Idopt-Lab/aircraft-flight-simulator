classdef SimpleMass < Mass
% SIMPLEMASS  Top-down mass model for conceptual design.
%
%   Specify empty mass, inertia, and CG directly via set_mass_properties().
%   Suitable for early design phases when component breakdown is unavailable.
%
%   Fuel tanks are handled separately and aggregated with empty properties
%   using the parallel axis theorem.
%
%   Example:
%     mass = SimpleMass();
%     mass.set_mass_properties(5000, 8000, 12000, 18000, 0, 0, 0, [2.5 0 0]);
%     mass.add_fuel_tank([3 0 0], 1500);
%     mass.set_fuel_mass(1200);
%     
%     euler = [0; deg2rad(5); 0];  % 5° pitch up
%     [F_g, M_g] = mass.get_gravity_force_moment(euler, [0; 0; 0]);
%
%   See also: Mass, ComponentMass, Aircraft

    properties
        % Empty aircraft mass [kg]
        m_empty = 0

        % Empty aircraft inertia tensor in body axes [kg-m^2]
        I_body_empty = zeros(3,3)

        % Empty aircraft CG location in body axes [m]
        cg_empty = [0 0 0]
    end

    methods

        function obj = SimpleMass()
        % SIMPLEMASS  Constructor. Initialises empty fuel tank array.

            obj.fuel_tanks = struct('location',{},'capacity',{},'current',{},'type',{});
        end

        function set_mass_properties(obj, m_empty, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, cg)
        % SET_MASS_PROPERTIES  Set empty-aircraft mass, inertia, and CG directly.
        %
        %   Inputs:
        %     m_empty - empty mass [kg]
        %     Ixx     - moment of inertia about x-axis [kg-m^2]
        %     Iyy     - moment of inertia about y-axis [kg-m^2]
        %     Izz     - moment of inertia about z-axis [kg-m^2]
        %     Ixy     - product of inertia xy [kg-m^2]
        %     Ixz     - product of inertia xz [kg-m^2]
        %     Iyz     - product of inertia yz [kg-m^2]
        %     cg      - CG position in body axes [x y z] [m]
        %
        %   Example:
        %     mass.set_mass_properties(5000, 8000, 12000, 18000, 0, 0, 0, [2.5 0 0]);

            obj.m_empty      = m_empty;
            obj.I_body_empty = [Ixx -Ixy -Ixz; -Ixy Iyy -Iyz; -Ixz -Iyz Izz];
            obj.cg_empty     = cg(:)';
            obj.update_current_properties();
        end

        % ── Required abstract implementations ────────────────────────────

        function m = get_total_mass(obj)
        % GET_TOTAL_MASS  Return total mass (empty + fuel) [kg].

            m = obj.m_empty + obj.fuel_mass;
        end

        function m = get_empty_mass(obj)
        % GET_EMPTY_MASS  Return empty mass [kg].

            m = obj.m_empty;
        end

        function cg = get_cg(obj)
        % GET_CG  Return current CG location in body axes [m].

            cg = obj.cg_current;
        end

        function I = get_inertia_matrix(obj)
        % GET_INERTIA_MATRIX  Return current inertia tensor [kg-m^2].

            I = obj.I_body_current;
            if rank(I) < 3
                error('SimpleMass:SingularInertia', ...
                    'Inertia matrix is singular. Call set_mass_properties first.');
            end
        end

    end

    methods (Access = protected)

        function update_current_properties(obj)
        % UPDATE_CURRENT_PROPERTIES  Recompute CG and inertia for current fuel state.

            m_total = obj.m_empty + obj.fuel_mass;
            if m_total < 1e-9
                obj.cg_current     = obj.cg_empty;
                obj.I_body_current = obj.I_body_empty;
                return;
            end

            % Fuel CG from tank locations
            if isempty(obj.fuel_tanks)
                fuel_cg = [0 0 0];
            else
                total_fuel  = 0;
                weighted_cg = [0 0 0];
                for i = 1:numel(obj.fuel_tanks)
                    m_t = obj.fuel_tanks(i).current;
                    total_fuel  = total_fuel  + m_t;
                    weighted_cg = weighted_cg + m_t * obj.fuel_tanks(i).location;
                end
                if total_fuel > 1e-9
                    fuel_cg = weighted_cg / total_fuel;
                else
                    fuel_cg = [0 0 0];
                end
            end

            % Combined CG
            obj.cg_current = (obj.m_empty*obj.cg_empty + obj.fuel_mass*fuel_cg) / m_total;

            % Empty inertia shifted to current CG
            dcg_empty = obj.cg_empty - obj.cg_current;
            I_empty_shifted = obj.parallel_axis_theorem(obj.I_body_empty, obj.m_empty, dcg_empty);

            % Fuel inertia contributions
            I_fuel_total = zeros(3,3);
            if ~isempty(obj.fuel_tanks)
                for i = 1:numel(obj.fuel_tanks)
                    m_t = obj.fuel_tanks(i).current;
                    if m_t < 1e-9, continue; end
                    dcg = obj.fuel_tanks(i).location - obj.cg_current;
                    if strcmpi(obj.fuel_tanks(i).type, 'distributed')
                        dx = dcg(1); dy = dcg(2); dz = dcg(3);
                        I_local = m_t * [dy^2+dz^2, 0, 0;
                                         0, dx^2+dz^2, 0;
                                         0, 0, dx^2+dy^2] / 12;
                    else
                        I_local = (2/5) * m_t * norm(dcg)^2 * eye(3);
                    end
                    I_fuel_total = I_fuel_total + obj.parallel_axis_theorem(I_local, m_t, dcg);
                end
            elseif obj.fuel_mass > 1e-9
                dcg = fuel_cg - obj.cg_current;
                if norm(dcg) > 1e-9
                    I_fuel_total = obj.parallel_axis_theorem((2/5)*obj.fuel_mass*norm(dcg)^2*eye(3), obj.fuel_mass, dcg);
                end
            end

            obj.I_body_current = I_empty_shifted + I_fuel_total;
        end

    end
end