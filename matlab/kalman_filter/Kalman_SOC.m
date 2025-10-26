clc
clear
close all 
%%Setup/Define Model Parameters
load("soc_ocv_data.mat"); 
SOC = SOC/100;
R_0 = 32e-3; %Internal Resistance Per Manufacturer
OCV = OCV + R_0*3; %Voltage/SOC curve provided gives terminal voltages at 1C. Approximate OCV = V_t + R_0*I
load("discharge_data.mat");
dt = zeros(size(data,1),1);

for i=1:size(data,1)-1
    dt(i+1) = data.time_s(i+1) - data.time_s  (i,1); 
end 

OCV_SOC = fit(SOC,OCV,'poly3'); %Fit Manufacture Provided Data to 3rd Order Poly
plot(SOC,OCV_SOC.p1*SOC.^3 + OCV_SOC.p2*SOC.^2 ...
        + OCV_SOC.p3.*SOC + OCV_SOC.p4);
hold on; 
plot(SOC,OCV,'o');


R_1 = 1.6e-3; %Diffusion and Charge Transfer Parameters/Per Literature as Placeholder for now
tau_1 = 3.68;
R_2 = 7.7e-3;
tau_2 = 84.34;
Q_nom = 3; %Capacity in mAmp-Hours

Q = 0.01; %Measurement Variance
R = diag([1000*Q;Q;0.1*Q]); %Process Covariance 

%%Extended Kalman Filter

x_predict = zeros(3,size(data,1));
x_update = zeros(3,size(data,1));
x = [100; 0; 0]; %Initial State Vector
SOC_0 = 1;

Sigma = eye(3); %Initial State Uncertainty
SOCs_Kalman = zeros(size(data,1),1); SOCs_CC = zeros(size(data,1),1);
for i=1:size(data,1)
    
    u = -data.I_mA(i)/1000; %Define positve current as Discharging
    A = diag([1;exp(-dt(i)/tau_1);exp(-dt(i)/tau_2)]);
    B = [-dt(i)/(Q_nom*3600);R_1*(1-exp(-dt(i)/tau_1));R_2*(1-exp(-dt(i)/tau_2))];
    x = A*x + B*u; %Prediction
    x_predict(:,i) = x;
    Sigma = A*Sigma*A.' + R; %Uncertianty
    C = [(3*OCV_SOC.p1*x(1)^2 + 2*OCV_SOC.p2*x(1) + OCV_SOC.p3), ... 
        -1, -1]; %Jacobian
    K = Sigma*C.'*inv(C*Sigma*C.'+ Q); %Kalman Gain 
    V_predict = (OCV_SOC.p1*x(1)^3 + OCV_SOC.p2*x(1)^2 ...
        + OCV_SOC.p3*x(1) + OCV_SOC.p4) - x(2) - x(3) - (u/1000)*R_0; 
    x = x + K*(data.Ecell_V(i) - V_predict); %Update
    x_update(:,i) = x;
    SOCs_Kalman(i) = x(1); 
    Sigma = (eye(3) - K*C)*Sigma; %Updated Uncertianty

end 



SOCs_CC = SOC_0 + cumtrapz(data.time_s,data.I_mA)/(Q_nom*3600*1000);
figure;

SOCs_Kalman = kalmanSOC(OCV_SOC,R_0,R_1,tau_1,R_2,tau_2,Q_nom,data); 

plot(data.time_s - data.time_s(1), SOCs_Kalman);
figure;
plot(data.time_s - data.time_s(1), SOCs_CC);
figure; 
plot(data.time_s - data.time_s(1),data.I_mA);

