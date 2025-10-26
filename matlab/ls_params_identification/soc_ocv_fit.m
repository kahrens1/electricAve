
clc;
clear;
close all

SOC = [1 0.95 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.25 0.2 0.15 0.1 0.05];
OCV = [4.17497 4.1042 4.05852 3.94657 3.86229 3.76835 3.66348 3.60236 ... 
    3.55024 3.51292 3.45824 3.39068 3.34436 3.23691];


% 3rd-Order Polynomial Fitting
degree = 8;
coefficients = polyfit(SOC, OCV, degree);

% Display coefficients
disp('Coefficients of the 3rd-order polynomial:');
disp(coefficients);

% Generate Fitted Curve
SOC_fit = linspace(0, 1, 100);
OCV_fit = polyval(coefficients, SOC_fit);

% Plotting
figure;
plot(SOC, OCV, 'o', 'MarkerSize', 8, 'DisplayName', 'Measured Data');
hold on;
plot(SOC_fit, OCV_fit, '-r', 'LineWidth', 2, 'DisplayName', 'Fitted Curve');
xlabel('SOC','FontSize', 16);
ylabel('OCV (V)','FontSize', 16);
%title('SOC-OCV Curve (8th-Order Polynomial)');
legend('show', 'Location', 'northwest','FontSize', 20);
grid on;

% Evaluate Fit Quality
OCV_predicted = polyval(coefficients, SOC);
mse = mean((OCV - OCV_predicted).^2);
disp(['Mean Squared Error: ', num2str(mse)]);