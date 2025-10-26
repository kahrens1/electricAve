clc;
clear;
close;

% Load data
load('HPPC_pulse_data_90.mat');

time = data.Time;
voltage = data.Voltage;
current = data.Current;

% Fixed OCV value
fixed_OCV = 4.05852;

% ECM model function
function [V_model] = ecm_model_fixed_ocv(params, current, time, fixed_OCV)
    % Decompose parameters
    R0 = params(1);
    R1 = params(2);
    C1 = params(3);
    R2 = params(4);
    C2 = params(5);

    % Time step
    dt = [0; diff(time)];

    % Initialize dynamic voltage components
    Vdl = 0; % Voltage across R1-C1
    Vdf = 0; % Voltage across R2-C2
    V_model = zeros(size(current));

    % Compute model voltage iteratively
    for k = 2:length(current)
        % Calculate exponential terms for RC branches
        a1 = exp(-dt(k) / (R1 * C1));
        b1 = R1 * (1 - exp(-dt(k) / (R1 * C1)));
        a2 = exp(-dt(k) / (R2 * C2));
        b2 = R2 * (1 - exp(-dt(k) / (R2 * C2)));

        % Update Vdl and Vdf based on current
        Vdl = a1 * Vdl + b1 * current(k - 1);
        Vdf = a2 * Vdf + b2 * current(k - 1);

        % Compute total model voltage
        V_model(k) = fixed_OCV + R0 * current(k) + Vdl + Vdf;
    end
end

% Error function for optimization
function error = ecm_error_fixed_ocv(params, current, voltage, time, fixed_OCV)
    V_model = ecm_model_fixed_ocv(params, current, time, fixed_OCV);
    error = sum((voltage - V_model).^2);
end

% Jacobian matrix computation
function J = compute_jacobian(params, current, time, delta, voltage, fixed_OCV)
    num_params = length(params);
    num_data = length(voltage);
    J = zeros(num_data, num_params);

    for j = 1:num_params
        params_forward = params;
        params_backward = params;
        params_forward(j) = params_forward(j) + delta;
        params_backward(j) = params_backward(j) - delta;

        V_forward = ecm_model_fixed_ocv(params_forward, current, time, fixed_OCV);
        V_backward = ecm_model_fixed_ocv(params_backward, current, time, fixed_OCV);

        % Numerical differentiation
        J(:, j) = (V_forward - V_backward) / (2 * delta);
    end
end

% Initial guesses and bounds
params0 = [0.02, 0.01, 100, 0.01, 200]; % Initial guesses: [R0, R1, C1, R2, C2]
lb = [0, 0, 1, 0, 1];  % Lower bounds
ub = [1, 10, 1e5, 10, 1e5]; % Upper bounds

% Optimization
options = optimset('Display', 'iter', 'TolFun', 1e-6, 'MaxIter', 500);
params_opt = fmincon(@(params) ecm_error_fixed_ocv(params, current, voltage, time, fixed_OCV), ...
                     params0, [], [], [], [], lb, ub, [], options);

% Extract optimized parameters
R0_opt = params_opt(1);
R1_opt = params_opt(2);
C1_opt = params_opt(3);
R2_opt = params_opt(4);
C2_opt = params_opt(5);

% Display results
fprintf('Optimal Parameters:\n');
fprintf('Fixed OCV = %.4f V\n', fixed_OCV);
fprintf('R0        = %.4f Ω\n', R0_opt);
fprintf('R1        = %.4f Ω\n', R1_opt);
fprintf('C1        = %.4f F\n', C1_opt);
fprintf('R2        = %.4f Ω\n', R2_opt);
fprintf('C2        = %.4f F\n', C2_opt);

% Calculate model voltage
V_model_opt = ecm_model_fixed_ocv(params_opt, current, time, fixed_OCV);

% Plot measured and modeled voltage
start_idx = 3;
end_idx = 600;

time_selected = time(start_idx:end_idx);
voltage_selected = voltage(start_idx:end_idx);
V_model_opt_selected = V_model_opt(start_idx:end_idx);

figure;
plot(time_selected, voltage_selected, 'b', 'LineWidth', 2, 'DisplayName', 'Measured Voltage');
hold on;
plot(time_selected, V_model_opt_selected, 'r--', 'LineWidth', 2, 'DisplayName', 'Model Voltage');
xlabel('Time (s)', 'FontSize', 20);
ylabel('Voltage (V)', 'FontSize', 20);
ylim([3.8 4.2]);
legend('FontSize', 20, 'Location', 'northeast');
grid on;
title('Measured Voltage vs. Estimated Voltage', 'FontSize', 20);

% Residual covariance matrix
residual = voltage - V_model_opt;
sigma_squared = var(residual);
delta = 1e-5; 
J = compute_jacobian(params_opt, current, time, delta, voltage, fixed_OCV);
param_covariance = sigma_squared * inv(J' * J);

% Display covariance matrix and standard deviations
disp('Parameter Covariance Matrix:');
disp(param_covariance);
disp('Parameter Standard Deviations:');
disp(sqrt(diag(param_covariance)));
