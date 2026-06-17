classdef AircraftConfigurator < handle

    properties
        aircraft
    end

    methods

        function obj = AircraftConfigurator(ac)
            obj.aircraft = ac;
        end

        % ================================================================
        % Generic component helper
        % ================================================================

        function comp = add_component(obj, varargin)

            p = inputParser;
            addParameter(p, 'name', '');
            addParameter(p, 'type', 'generic');
            addParameter(p, 'mass', 0);
            addParameter(p, 'position', [0 0 0]);
            addParameter(p, 'cg_local', [0 0 0]);
            addParameter(p, 'inertia', zeros(3,3));
            addParameter(p, 'parent_frame', obj.aircraft.body_frame_name);
            addParameter(p, 'frame_name', '');
            addParameter(p, 'dcm_fn', @(x) eye(3));
            addParameter(p, 'parent_component', '');
            parse(p, varargin{:});
            s = p.Results;

            if strlength(string(s.frame_name)) > 0
                frame_name = char(s.frame_name);
            else
                frame_name = matlab.lang.makeValidName(char(string(s.name) + "_frame"));
            end

            if ~obj.aircraft.has_frame(frame_name)
                obj.aircraft.add_frame( ...
                    frame_name, ...
                    s.parent_frame, ...
                    s.position(:), ...
                    s.dcm_fn);
            else
                obj.aircraft.update_frame_position(frame_name, s.position(:));
                obj.aircraft.update_frame_orientation(frame_name, s.dcm_fn);
            end

            comp_frame = obj.aircraft.get_frame(frame_name);

            comp = Component( ...
                s.name, ...
                s.mass, ...
                s.cg_local, ...
                s.inertia, ...
                comp_frame, ...
                s.type);

            if strlength(string(s.parent_component)) > 0
                obj.aircraft.add_component(comp, s.parent_component);
            else
                obj.aircraft.add_component(comp);
            end
        end

        function comp = add_mass_component(obj, varargin)
            % Backward-compatible wrapper.
            comp = obj.add_component(varargin{:});
        end

        % ================================================================
        % Control surfaces
        % ================================================================

        function add_control_surface(obj, varargin)

            p = inputParser;
            addParameter(p, 'name', '');
            addParameter(p, 'surface_type', '');
            addParameter(p, 'classification', '');
            addParameter(p, 'axis', [0 0 0]);
            addParameter(p, 'max_deflection', 0);
            addParameter(p, 'min_deflection', 0);
            addParameter(p, 'dCl', 0);
            addParameter(p, 'dCm', 0);
            addParameter(p, 'dCn', 0);
            parse(p, varargin{:});
            s = p.Results;

            cs = ControlSurface( ...
                s.name, ...
                s.surface_type, ...
                s.classification, ...
                s.axis, ...
                s.max_deflection, ...
                s.min_deflection, ...
                s.dCl, ...
                s.dCm, ...
                s.dCn);

            obj.aircraft.add_control_surface(cs);
        end

        % ================================================================
        % Propulsion
        % ================================================================

        function [pe, comp] = add_propulsive_element(obj, varargin)

            p = inputParser;
            addParameter(p, 'name', '');
            addParameter(p, 'element_type', '');
            addParameter(p, 'max_output', 0);
            addParameter(p, 'position', [0 0 0]);
            addParameter(p, 'direction', [1 0 0]);
            addParameter(p, 'mount_euler', [0 0 0]);
            addParameter(p, 'fuel_rate', 0);
            addParameter(p, 'thrust_model', []);

            % Component/mass options
            addParameter(p, 'mass', 0);
            addParameter(p, 'inertia', zeros(3,3));
            addParameter(p, 'cg_local', [0 0 0]);
            addParameter(p, 'parent_component', '');
            addParameter(p, 'parent_frame', obj.aircraft.body_frame_name);

            % Load options
            addParameter(p, 'add_load_source', true);
            addParameter(p, 'frame_name', '');
            parse(p, varargin{:});
            s = p.Results;

            dir = s.direction(:);

            if norm(dir) < 1e-9
                dir = [1;0;0];
            else
                dir = dir / norm(dir);
            end

            if strlength(string(s.frame_name)) > 0
                frame_name = char(s.frame_name);
            else
                frame_name = matlab.lang.makeValidName(char(string(s.name) + "_thrust_frame"));
            end

            if ~obj.aircraft.has_frame(frame_name)
                obj.aircraft.add_frame( ...
                    frame_name, ...
                    s.parent_frame, ...
                    s.position(:), ...
                    @(x) obj.mount_dcm_from_euler(s.mount_euler));
            else
                obj.aircraft.update_frame_position(frame_name, s.position(:));
                obj.aircraft.update_frame_orientation(frame_name, ...
                    @(x) obj.mount_dcm_from_euler(s.mount_euler));
            end

            thrust_frame = obj.aircraft.get_frame(frame_name);

            et = lower(strtrim(char(s.element_type)));

            switch et
                case {'turbojet','turbofan','turbofan_afterburning'}
                    pe = TurbofanPropulsion( ...
                        s.name, ...
                        thrust_frame, ...
                        dir, ...
                        s.max_output, ...
                        s.fuel_rate, ...
                        s.thrust_model);

                case {'propeller','prop'}
                    pe = PropellerPropulsion( ...
                        s.name, ...
                        thrust_frame, ...
                        dir, ...
                        s.max_output, ...
                        s.fuel_rate, ...
                        s.thrust_model);

                case {'rotor','rotor_variable','multirotor','rotor_blockset'}
                    pe = RotorPropulsion( ...
                        s.name, ...
                        thrust_frame, ...
                        dir, ...
                        s.thrust_model);

                case {'electric'}
                    pe = ElectricPropulsion( ...
                        s.name, ...
                        thrust_frame, ...
                        dir, ...
                        s.max_output, ...
                        s.thrust_model);

                case {'directforce','direct_force','bodyforce'}
                    pe = DirectForcePropulsiveElement( ...
                        s.name, ...
                        thrust_frame, ...
                        dir, ...
                        s.thrust_model);

                otherwise
                    pe = CustomPropulsion( ...
                        s.name, ...
                        thrust_frame, ...
                        dir, ...
                        s.thrust_model);
            end

            obj.aircraft.add_propulsive_element(pe);

            comp = Component( ...
                s.name, ...
                s.mass, ...
                s.cg_local, ...
                s.inertia, ...
                thrust_frame, ...
                "engine");

            if s.add_load_source
                pls = PropulsionLoadSolver(pe, thrust_frame);
                comp.add_load_source(pls);
            end

            if strlength(string(s.parent_component)) > 0
                obj.aircraft.add_component(comp, s.parent_component);
            else
                obj.aircraft.add_component(comp);
            end
        end

        % ================================================================
        % Aerodynamics
        % ================================================================

        function solver = add_aero_solver(obj, aero_model, varargin)

            p = inputParser;
            addParameter(p, 'name', 'aircraft_aero');
            addParameter(p, 'position', obj.aircraft.geometry.ref_point);
            addParameter(p, 'dcm_fn', @(x) eye(3));
            addParameter(p, 'geom', obj.aircraft.geometry);
            addParameter(p, 'component_name', '');  % empty = root aircraft component
            addParameter(p, 'parent_frame', obj.aircraft.body_frame_name);
            addParameter(p, 'frame_name', '');
            parse(p, varargin{:});
            s = p.Results;

            if strlength(string(s.frame_name)) > 0
                frame_name = char(s.frame_name);
            else
                frame_name = matlab.lang.makeValidName(char(string(s.name) + "_frame"));
            end

            if ~obj.aircraft.has_frame(frame_name)
                obj.aircraft.add_frame( ...
                    frame_name, ...
                    s.parent_frame, ...
                    s.position(:), ...
                    s.dcm_fn);
            else
                obj.aircraft.update_frame_position(frame_name, s.position(:));
                obj.aircraft.update_frame_orientation(frame_name, s.dcm_fn);
            end

            aero_frame = obj.aircraft.get_frame(frame_name);

            solver = AeroLoadSolver( ...
                aero_model, ...
                s.geom, ...
                obj.aircraft, ...
                aero_frame);

            if strlength(string(s.component_name)) > 0
                obj.aircraft.add_load_source(solver, s.component_name);
            else
                obj.aircraft.add_load_source(solver);
            end
        end

        % ================================================================
        % Printing
        % ================================================================

        function print_configuration(obj)

            ac = obj.aircraft;

            fprintf('\n=== AIRCRAFT CONFIGURATION ===\n');

            fprintf('\nFrames: %d\n', numel(ac.list_frames()));
            frame_names = ac.list_frames();

            for i = 1:numel(frame_names)
                fprintf('  %d) %s\n', i, char(frame_names(i)));
            end

            fprintf('\nComponent Tree:\n');

            if ismethod(ac,'print_component_tree')
                ac.print_component_tree();
            else
                for k = 1:numel(ac.components)
                    c = ac.components{k};
                    fprintf('  %d) %-15s mass=%.3f kg\n', k, char(c.name), c.mass);
                end
            end

            fprintf('\nControl Surfaces: %d\n', numel(ac.control_surfaces));

            for i = 1:numel(ac.control_surfaces)
                cs = ac.control_surfaces(i);

                fprintf('  %d) %-12s type=%-10s class=%-10s axis=[%.0f %.0f %.0f]\n', ...
                    i, ...
                    char(string(cs.name)), ...
                    char(string(cs.surface_type)), ...
                    char(string(cs.classification)), ...
                    cs.axis(1), ...
                    cs.axis(2), ...
                    cs.axis(3));
            end

            fprintf('\nPropulsive Elements: %d\n', numel(ac.propulsive_elements));

            for k = 1:numel(ac.propulsive_elements)

                pe = ac.propulsive_elements{k};

                pe_name  = "";
                pe_class = class(pe);
                pe_axis  = [1;0;0];
                pe_frame = "none";

                if isprop(pe, 'name')
                    pe_name = string(pe.name);
                end

                if isprop(pe, 'local_thrust_axis')
                    pe_axis = pe.local_thrust_axis(:);
                end

                if isprop(pe, 'frame') && ~isempty(pe.frame)
                    pe_frame = string(pe.frame.name);
                end

                fprintf('  %d) %-15s class=%-25s frame=%-20s axis=[%.3f %.3f %.3f]\n', ...
                    k, ...
                    char(pe_name), ...
                    pe_class, ...
                    char(pe_frame), ...
                    pe_axis(1), pe_axis(2), pe_axis(3));
            end

            fprintf('==============================\n\n');
        end

    end

    methods (Static)

        function C = mount_dcm_from_euler(eul)

            eul = eul(:);

            if numel(eul) < 3
                eul(end+1:3,1) = 0;
            end

            roll  = eul(1);
            pitch = eul(2);
            yaw   = eul(3);

            cr = cos(roll);  sr = sin(roll);
            cp = cos(pitch); sp = sin(pitch);
            cy = cos(yaw);   sy = sin(yaw);

            Rx = [1 0 0;
                  0 cr -sr;
                  0 sr  cr];

            Ry = [ cp 0 sp;
                    0 1 0;
                  -sp 0 cp];

            Rz = [cy -sy 0;
                  sy  cy 0;
                   0   0 1];

            C = Rz * Ry * Rx;
        end

    end
end