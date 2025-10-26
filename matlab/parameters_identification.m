clc;
clear;
close all

% Load data from .mat file
load('HPPC_pulse_data_30.mat'); 

% Load the data
time = data.Time; 
voltage = data.Voltage; 
current = data.Current; 


% ECM model function
function [V_model] = ecm_model(params, current, time)
    % Decompose parameters
    OCV = params(1); 
    R0 = params(2);
    R1 = params(3); 
    C1 = params(4); 
    R2 = params(5); 
    C2 = params(6);
    
    % Time step
    dt = [0; diff(time)]; 
    
    % Initialize dynamic voltage components
    Vdl = 0; % Initial Vdl (voltage across R1-C1 branch)
    Vdf = 0; % Initial Vdf (voltage across R2-C2 branch)
    V_model = zeros(size(current)); % Preallocate model voltage array
    
    % Compute model voltage iteratively
    for k = 2:length(current)
        % Calculate exponential terms for RC branches
        a1 = exp(-dt(k) / (R1 * C1));
        b1 = R1 * (1 - exp(-dt(k) / (R1 * C1)));
        a2 = exp(-dt(k) / (R2 * C2));
        b2 = R2 * (1 - exp(-dt(k) / (R2 * C2)));
        
        % Update Vdl and Vdf based on current
        Vdl = a1 * Vdl + b1 * current(k-1);
        Vdf = a2 * Vdf + b2 * current(k-1);
        
        % Compute the total model voltage
        V_model(k) = OCV + R0 * current(k) + Vdl + Vdf;
    end
end

% Error function for parameter optimization
function error = ecm_error(params, current, voltage, time)
    % Compute the model voltage with current parameters
    V_model = ecm_model(params, current, time);
    % Compute the sum of squared errors between model and measured voltage
    error = sum((voltage - V_model).^2); 
end

% Initial parameter estimates
params0 = [3.6, 0.2, 0.05, 100, 0.05, 200]; % Initial guesses: [OCV, R0, R1, C1, R2, C2]

% Optimization options
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'MaxIter', 500 );

% Perform optimization to minimize error
params_opt = fmincon(@(params) ecm_error(params, current, voltage, time), ...
                     params0, [], [], [], [], ...
                     [3, 0, 0, 0, 0, 0], ...  % Lower bounds for parameters
                     [4.5, Inf, Inf, Inf, Inf, Inf], ... % Upper bounds for parameters
                     [], options);

% Extract optimized parameters
OCV_opt = params_opt(1);
R0_opt = params_opt(2);
R1_opt = params_opt(3);
C1_opt = params_opt(4);
R2_opt = params_opt(5);
C2_opt = params_opt(6);

% Display optimized parameters
fprintf('Optimal Parameters:\n');
fprintf('OCV   = %.4f V\n', OCV_opt);
fprintf('R0    = %.4f Ω\n', R0_opt);
fprintf('R1    = %.4f Ω\n', R1_opt);
fprintf('C1    = %.4f F\n', C1_opt);
fprintf('R2    = %.4f Ω\n', R2_opt);
fprintf('C2    = %.4f F\n', C2_opt);

% Compute the optimized model voltage
V_model_opt = ecm_model(params_opt, current, time);

% Plot measured and modeled voltage
figure;
plot(time, voltage, 'b', 'DisplayName', 'Measured Voltage');
hold on;
plot(time, V_model_opt, 'r--', 'DisplayName', 'Model Voltage');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend;
grid on;
title('ECM Fitting Results with OCV Optimization');