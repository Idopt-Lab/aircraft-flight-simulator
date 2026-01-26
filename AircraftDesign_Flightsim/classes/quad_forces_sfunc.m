function quad_forces_sfunc(block)
setup(block);
end

function setup(block)
block.NumInputPorts  = 2;
block.NumOutputPorts = 2;

block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions        = 12;
block.InputPort(1).DirectFeedthrough = true;

block.InputPort(2).Dimensions        = 4;
block.InputPort(2).DirectFeedthrough = true;

block.OutputPort(1).Dimensions = 3;
block.OutputPort(2).Dimensions = 3;

block.SampleTimes        = [0 0];
block.SimStateCompliance = 'DefaultSimState';

block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Outputs', @Outputs);
end

function DoPostPropSetup(block)
block.NumDworks = 0;
end

function InitializeConditions(block)
end

function Outputs(block)
persistent ac n_cs n_pe

if isempty(ac)
    try
        ac = evalin('base','AIRCRAFT_OBJ');
    catch
        try
            ac = evalin('base','quad');
        catch
            try
                ac = evalin('base','ac');
            catch
                error('Cannot find aircraft object in base workspace');
            end
        end
    end
    
    n_cs = numel(ac.control_surfaces);
    n_pe = numel(ac.propulsive_elements);
end

x = block.InputPort(1).Data(:);
u = block.InputPort(2).Data(:);

ac.state.set_full_state(x);

for i = 1:n_cs
    if i <= numel(u)
        ac.control_surfaces(i).set_deflection(u(i));
    else
        ac.control_surfaces(i).set_deflection(0);
    end
end

for k = 1:n_pe
    j = n_cs + k;
    if j <= numel(u)
        throttle_cmd = max(0, min(1, u(j)));
        ac.propulsive_elements{k}.set_throttle(throttle_cmd);
    else
        ac.propulsive_elements{k}.set_throttle(0);
    end
end

ac.sync_control_vector_from_components();

[F_ext, M_ext, ~] = ac.calculate_external_forces_moments();

F_ext = F_ext(:);
M_ext = M_ext(:);

if ~all(isfinite(F_ext))
    F_ext = zeros(3,1);
end

if ~all(isfinite(M_ext))
    M_ext = zeros(3,1);
end

block.OutputPort(1).Data = F_ext;
block.OutputPort(2).Data = M_ext;
end