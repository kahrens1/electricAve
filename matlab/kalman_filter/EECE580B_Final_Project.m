clc
clear 
close all 

%%Setup/Define Model Parameters
load("soc_ocv_data.mat"); 
SOC = SOC/100;
R_0 = 32e-3; %Internal Resistance Per Manufacturer
OCV = OCV + R_0*3; %Voltage/SOC curve provided gives terminal voltages at 1C. Approximate OCV = V_t + R_0*I
load("discharge_data.mat");

OCV_SOC = fit(SOC,OCV,'poly3'); %Fit Manufacture Provided Data to 3rd Order Poly
OCV_SOC_fit = OCV_SOC.p1*SOC.^3 + OCV_SOC.p2*SOC.^2 ...
        + OCV_SOC.p3.*SOC + OCV_SOC.p4;
subplot(2,1,1);
plot(SOC,OCV_SOC.p1*SOC.^3 + OCV_SOC.p2*SOC.^2 ...
        + OCV_SOC.p3.*SOC + OCV_SOC.p4);
xlabel('SOC (%)');
ylabel("Open Circuit Voltage (V)");

hold on;
plot(SOC,OCV,'o');
subplot(2,1,2);
plot(SOC, abs(OCV - OCV_SOC_fit));



R_1 = 1.6e-3; %Diffusion and Charge Transfer Parameters/Per Literature as Placeholder for now
tau_1 = 3.68;
R_2 = 7.7e-3;
tau_2 = 84.34;
Q_nom = 3; %Capacity in Amp-Hours
SOC_0 = 1;

SOCs_Kalman_Dischage = kalmanSOC(OCV_SOC,R_0,R_1,tau_1,R_2,tau_2,Q_nom,data); 
figure;
plot(data.time_s - data.time_s(1), SOCs_Kalman_Dischage*100);
SOCs_CC = 100*(SOC_0 + cumtrapz(data.time_s,data.I_mA)/(Q_nom*3600*1000));
hold on
plot(data.time_s -data.time_s(1),SOCs_CC);
legend('Kalman Filter','Coloumb Counting');
xlabel('Time (s)');
ylabel('State of Charge (%)');

load('charge_data.mat'); 
SOCs_Kalman_Charge = 100*(kalmanSOC(OCV_SOC,R_0,R_1,tau_1,R_2,tau_2,Q_nom,data)); 
figure;
plot(data.time_s - data.time_s(1),SOCs_Kalman_Charge);
hold on 
SOC_0 = 0.25; %Approximate (initial SOC not stated in data)
SOCs_CC = 100*(SOC_0 + cumtrapz(data.time_s,data.I_mA)/(Q_nom*3600*1000));
plot(data.time_s - data.time_s(1),SOCs_CC);
legend('Kalman Filter','Coloumb Counting');
xlabel('Time (s)');
ylabel('State of Charge (%)');


load('charge_discharge_data.mat'); 
SOCs_Kalman_Charge_Discharge = 100*(kalmanSOC(OCV_SOC,R_0,R_1,tau_1,R_2,tau_2,Q_nom,data)); 
figure;
plot(data.time_s - data.time_s(1),SOCs_Kalman_Charge_Discharge);
hold on 
SOC_0 = 1; %Approximate (initial SOC not stated in data)
SOCs_CC = 100*(SOC_0 + cumtrapz(data.time_s,data.I_mA)/(Q_nom*3600*1000));
plot(data.time_s - data.time_s(1),SOCs_CC);
legend('Kalman Filter','Coloumb Counting');
xlabel('Time (s)');
ylabel('State of Charge (%)');