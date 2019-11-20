armParamHW8;  % load parameters

% instantiate arm, controller, and reference input classes 
% Instantiate Dynamics class
addpath('../hw3'); arm = armDynamics(P);  
controller = armController(P);  
addpath('../hw2'); reference = signalGenerator(30*pi/180, 0.05);  
addpath('../hw2'); disturbance = signalGenerator(0.0, 0.0);

% instantiate the data plots and animation
addpath('../hw2'); dataPlot = plotData(P);
addpath('../hw2'); animation = armAnimation(P);

% main simulation loop
t = P.t_start;  % time starts at t_start
y = arm.h();  % output at time t_start
while t < P.t_end  
    t_next_plot = t + P.t_plot; % only update plots at rate t_plot
    while t < t_next_plot 
        r = reference.square(t);  % compute reference
        d = disturbance.step(t);  % compute disturbance
        x = arm.state;  % use state feedback instead of output
        u = controller.update(r, x);  % Calculate the control
        y = arm.update(u+d);  % Propagate the dynamics
        t = t + P.Ts; % advance time by Ts
    end
    % update animation and data plots
    animation.update(arm.state);
    dataPlot.update(t, r, arm.state, u);
end


