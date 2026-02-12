function quad_forces_sfunc(block)
setup(block);
end

function setup(block)
block.NumInputPorts  = 2;
block.NumOutputPorts = 3;
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions        = 12;
block.InputPort(1).DirectFeedthrough = true;
block.InputPort(2).Dimensions        = 1;
block.InputPort(2).DirectFeedthrough = false;

block.OutputPort(1).Dimensions = 3;
block.OutputPort(2).Dimensions = 3;
block.OutputPort(3).Dimensions = 1;

block.SampleTimes = [0.01 0];
block.SimStateCompliance = 'DefaultSimState';

block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Outputs', @Outputs);
end

function DoPostPropSetup(block)
block.NumDworks = 2;

block.Dwork(1).Name = 'h_integral';
block.Dwork(1).Dimensions = 1;
block.Dwork(1).DatatypeID = 0;
block.Dwork(1).Complexity = 'Real';

block.Dwork(2).Name = 'last_t';
block.Dwork(2).Dimensions = 1;
block.Dwork(2).DatatypeID = 0;
block.Dwork(2).Complexity = 'Real';
end

function InitializeConditions(block)
block.Dwork(1).Data = 0;
block.Dwork(2).Data = 0;
end

function Outputs(block)
persistent ac n_cs n_pe auto_mode

if isempty(ac)
    try
        ac = evalin('base','quad');
    catch
        try
            ac = evalin('base','ac');
        catch
            error('Cannot find aircraft object in base workspace');
        end
    end
    n_cs = numel(ac.control_surfaces);
    n_pe = numel(ac.propulsive_elements);
    auto_mode = 1;
end

x = block.InputPort(1).Data(:);
manual_throttle = block.InputPort(2).Data(1);

ac.state.set_full_state(x);

if evalin('base','exist(''auto_mode'',''var'')')
    auto_mode = evalin('base','auto_mode');
end

if auto_mode
    Kp = 2.5;
    Ki = 0.1;
    Kd = 3.0;
    u_hover = 0.613;
    
    if evalin('base','exist(''Kp_h'',''var'')')
        Kp = evalin('base','Kp_h');
    end
    if evalin('base','exist(''Ki_h'',''var'')')
        Ki = evalin('base','Ki_h');
    end
    if evalin('base','exist(''Kd_h'',''var'')')
        Kd = evalin('base','Kd_h');
    end
    if evalin('base','exist(''u_hover'',''var'')')
        u_hover = evalin('base','u_hover');
    end
    
    t = block.CurrentTime;
    
    if evalin('base','exist(''h_cmd_data'',''var'')')
        h_cmd_data = evalin('base','h_cmd_data');
        h_cmd = interp1(h_cmd_data(:,1), h_cmd_data(:,2), t, 'linear', 'extrap');
    else
        h_cmd = 10.0;
        if evalin('base','exist(''h_cmd'',''var'')')
            h_cmd = evalin('base','h_cmd');
        end
    end
    
    last_t = block.Dwork(2).Data;
    dt = t - last_t;
    if dt <= 0 || ~isfinite(dt)
        dt = 0.01;
    end
    block.Dwork(2).Data = t;
    
    h_current = -x(3);
    h_error = h_cmd - h_current;
    
    h_integral = block.Dwork(1).Data;
    h_integral = h_integral + h_error * dt;
    h_integral = max(min(h_integral, 5), -5);
    block.Dwork(1).Data = h_integral;
    
    hdot = -x(6);
    
    u_cmd = u_hover + Kp * h_error + Ki * h_integral - Kd * hdot;
    throttle_cmd = max(0.1, min(0.95, u_cmd));
else
    throttle_cmd = max(0, min(1, manual_throttle));
end

u_total = zeros(n_cs + n_pe, 1);

if n_pe > 0
    u_total(n_cs + 1) = throttle_cmd;
    ac.propulsive_elements{1}.set_throttle(throttle_cmd);
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
block.OutputPort(3).Data = throttle_cmd;
end