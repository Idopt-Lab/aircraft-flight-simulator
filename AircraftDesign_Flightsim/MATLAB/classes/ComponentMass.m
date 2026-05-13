classdef ComponentMass < Mass
    % COMPONENTMASS  Bottom-up aircraft mass model using discrete components.
    %
    %   Builds aircraft mass properties by aggregating individual structural,
    %   payload, systems, and fuel components.
    %
    %   Total CG and inertia are computed using mass-weighted averaging and
    %   the parallel axis theorem.
    %
    %   Modeling assumptions:
    %     1. All component positions are defined in aircraft body axes [m].
    %     2. Component inertia tensors are specified about each component CG.
    %     3. Total aircraft inertia is assembled about the current aircraft CG.
    %     4. Fuel tanks are represented using simplified lumped or approximate
    %        distributed-mass models.
    %     5. Fuel slosh and dynamic fuel migration are not modeled.
    %
    %   References:
    %     Stevens, B. L., Lewis, F. L., and Johnson, E. N.,
    %     Aircraft Control and Simulation, 3rd ed., Wiley, 2015.
    %
    %     Etkin, B. and Reid, L. D.,
    %     Dynamics of Flight: Stability and Control, 3rd ed., Wiley, 1996.
    %
    %   See also:
    %     Mass, SimpleMass, Aircraft

    properties
        % Struct array of mass components: name, mass, position, I_local, type
        components = []
    end

    methods

        function obj = ComponentMass()
            % COMPONENTMASS  Constructor. Initialises empty component and fuel arrays.

            obj.components = struct('name',{},'mass',{},'position',{},'I_local',{},'type',{});
            obj.fuel_tanks = struct('location',{},'capacity',{},'current',{},'type',{});
        end

        function add_component(obj, name, mass, position, I_local, comp_type)
            % ADD_COMPONENT  Add a mass component with position and inertia.
            %
            %   Inputs:
            %     name     - component identifier string
            %     mass     - component mass [kg]
            %     position - 1x3 or 3x1 body-frame position [m]
            %     I_local  - 3x3 inertia tensor about component CG [kg-m^2]
            %                (optional, default: point mass)
            %     comp_type- component type string (optional, e.g. 'structure',
            %                'avionics', 'payload')
            %
            %   Example:
            %     mass.add_component('fuselage', 500, [2 0 0], eye(3)*10, 'structure');

            if nargin < 5 || isempty(I_local)
                I_local = zeros(3,3);
            end
            if nargin < 6 || isempty(comp_type)
                comp_type = 'generic';
            end

            n = numel(obj.components);
            obj.components(n+1).name     = char(name);
            obj.components(n+1).mass     = mass;
            obj.components(n+1).position = position(:)';
            obj.components(n+1).I_local  = I_local;
            obj.components(n+1).type     = char(comp_type);

            obj.update_current_properties();
        end

        function remove_component(obj, name)
            % REMOVE_COMPONENT  Remove a component by name.
            %
            %   Input:
            %     name - component name string

            if isempty(obj.components), return; end
            idx = find(strcmp({obj.components.name}, char(name)), 1);
            if ~isempty(idx)
                obj.components(idx) = [];
                obj.update_current_properties();
            end
        end

        function list_components(obj)
            % LIST_COMPONENTS  Print a summary of all mass components.

            if ~isempty(obj.components)
                fprintf('\n=== MASS COMPONENTS ===\n');
                fprintf('%-15s %10s %18s %12s\n', 'Name', 'Mass [kg]', 'Position [m]', 'Type');
                fprintf('%s\n', repmat('-',60,1));
                for i = 1:numel(obj.components)
                    c = obj.components(i);
                    fprintf('%-15s %10.2f  [%5.2f %5.2f %5.2f]  %-12s\n', ...
                        c.name, c.mass, c.position(1), c.position(2), c.position(3), c.type);
                end
                fprintf('Total empty mass: %.2f kg\n', obj.get_empty_mass());
                fprintf('Current total mass: %.2f kg\n', obj.get_total_mass());
                fprintf('Current CG: [%.3f  %.3f  %.3f] m\n', obj.cg_current(1), obj.cg_current(2), obj.cg_current(3));
            else
                fprintf('No components defined.\n');
            end
        end

        % ── Required abstract implementations ────────────────────────────

        function m = get_total_mass(obj)
            % GET_TOTAL_MASS  Return total mass (empty + fuel) [kg].

            m = obj.get_empty_mass() + obj.fuel_mass;
        end

        function m = get_empty_mass(obj)
            % GET_EMPTY_MASS  Return sum of component masses [kg].

            m = 0;
            for i = 1:numel(obj.components)
                m = m + obj.components(i).mass;
            end
        end

        function cg = get_cg(obj)
            % GET_CG  Return current CG location in body axes [m].

            cg = obj.cg_current;
        end

        function I = get_inertia_matrix(obj)
            % GET_INERTIA_MATRIX  Return current inertia tensor [kg-m^2].

            I = obj.I_body_current;
            if rcond(I) < 1e-12
                error('ComponentMass:SingularInertia', ...
                    'Inertia matrix is singular. Add components first.');
            end
        end

    end

    methods (Access = protected)

        function update_current_properties(obj)
            % UPDATE_CURRENT_PROPERTIES  Aggregate components and fuel into total mass/CG/inertia.

            % Step 1: Aggregate empty components
            m_empty_total = 0;
            cg_weighted   = [0 0 0];
            for i = 1:numel(obj.components)
                m_i = obj.components(i).mass;
                m_empty_total = m_empty_total + m_i;
                cg_weighted   = cg_weighted + m_i * obj.components(i).position;
            end

            cg_empty = [0 0 0];
            if m_empty_total > 1e-9
                cg_empty = cg_weighted / m_empty_total;
            end

            % Aggregate empty inertia about empty CG
            I_empty_agg = zeros(3,3);
            for i = 1:numel(obj.components)
                m_i   = obj.components(i).mass;
                pos_i = obj.components(i).position;
                I_i   = obj.components(i).I_local;
                dcg   = pos_i - cg_empty;
                I_empty_agg = I_empty_agg + obj.parallel_axis_theorem(I_i, m_i, dcg);
            end

            % Step 2: Add fuel contribution
            m_total = m_empty_total + obj.fuel_mass;
            if m_total < 1e-9
                obj.cg_current     = cg_empty;
                obj.I_body_current = I_empty_agg;
                return;
            end

            % Fuel CG
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
            obj.cg_current = (m_empty_total*cg_empty + obj.fuel_mass*fuel_cg) / m_total;

            % Empty inertia shifted to current CG
            dcg_empty = cg_empty - obj.cg_current;
            I_empty_shifted = obj.parallel_axis_theorem(I_empty_agg, m_empty_total, dcg_empty);

            % Fuel inertia contributions
            I_fuel_total = zeros(3,3);
            if ~isempty(obj.fuel_tanks)
                for i = 1:numel(obj.fuel_tanks)
                    m_t = obj.fuel_tanks(i).current;
                    if m_t < 1e-9, continue; end
                    dcg = obj.fuel_tanks(i).location - obj.cg_current;
                    if strcmpi(obj.fuel_tanks(i).type, 'distributed')
                        % Simplified distributed-fuel approximation:
                        % intrinsic tank inertia neglected; tank placement handled using
                        % parallel-axis shift relative to aircraft CG.
                        I_local = zeros(3,3);
                    else
                        % Lumped/spherical approximation
                        I_local = zeros(3,3);
                    end
                    I_fuel_total = I_fuel_total + obj.parallel_axis_theorem(I_local, m_t, dcg);
                end
            elseif obj.fuel_mass > 1e-9
                % Fuel mass exists but no tank geometry is defined.
                % Treat fuel as a lumped point mass at fuel_cg.
                dcg = fuel_cg - obj.cg_current;
                I_fuel_total = obj.parallel_axis_theorem(zeros(3,3), obj.fuel_mass, dcg);
            end

            obj.I_body_current = I_empty_shifted + I_fuel_total;
        end

    end
end