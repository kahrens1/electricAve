function SOCs = kalmanSOC(OCV_SOC,R_0s,R_1s,tau_1s,R_2s,tau_2s,Q_nom,data)
      
assert(isa(OCV_SOC, 'cfit'), 'Input must be 3rd Order SOC/OCV Curve Approximation');
  
%%Extended Kalman Filter

%x_predict = zeros(3,size(data,1));
%x_update = zeros(3,size(data,1));
x = [0.5; 0.1; 0.1]; %Initial State Vector
%SOC_0 = 1;
Sigma = eye(3); %Initial State Uncertainty
SOCs = zeros(size(data,1),1);
Q = 0.01; %Measurement Variance
R = diag([1000*Q;Q;0.1*Q]); %Process Covariance 

dt = zeros(size(data,1),1);

for i=1:size(data,1)-1
    dt(i+1) = data.time_s(i+1) - data.time_s  (i,1); 
end 


for i=1:size(data,1)
    %Interpolate within SOC Ranges to Find Parameters 
    dSOC = 0.2;
   
    if  x(1) >= 0.9
        R_0 = R_0s(1);
        R_1 = R_1s(1);
        tau_1 = tau_1s(1);
        R_2 = R_2s(1);
        tau_2 = tau_2s(1);
    elseif (x(1) < 0.9) && (x(1) >= 0.7)
        j = 1;
        R_0 = ((R_0s(j) - R_0s(j+1)) / dSOC) * (x(1) - 0.9) + R_0s(j);
        R_1 = ((R_1s(j) - R_1s(j+1)) / dSOC) * (x(1) - 0.9) + R_1s(j);
        tau_1 = ((tau_1s(j) - tau_1s(j+1)) / dSOC) * (x(1) - 0.9) +tau_1s(j);
        R_2 = ((R_2s(j) - R_2s(j+1)) / dSOC) * (x(1) - 0.9) + R_2s(j);
        tau_2 = ((tau_2s(j) - tau_2s(j+1)) / dSOC) * (x(1) - 0.9) +tau_2s(j);
    elseif (x(1) < 0.7) && (x(1) >= 0.5)
        j = 2;
        R_0 = ((R_0s(j) - R_0s(j+1)) / dSOC) * (x(1) - 0.7) + R_0s(j);
        R_1 = ((R_1s(j) - R_1s(j+1)) / dSOC) * (x(1) - 0.7) + R_1s(j);
        tau_1 = ((tau_1s(j) - tau_1s(j+1)) / dSOC) * (x(1) - 0.7) +tau_1s(j);
        R_2 = ((R_2s(j) - R_2s(j+1)) / dSOC) * (x(1) - 0.7) + R_2s(j);
        tau_2 = ((tau_2s(j) - tau_2s(j+1)) / dSOC) * (x(1) - 0.7) +tau_2s(j);
    elseif (x(1) < 0.5) && (x(1) >= 0.3)
        j = 3;
        R_0 = ((R_0s(j) - R_0s(j+1)) / dSOC) * (x(1) - 0.5) + R_0s(j);
        R_1 = ((R_1s(j) - R_1s(j+1)) / dSOC) * (x(1) - 0.5) + R_1s(j);
        tau_1 = ((tau_1s(j) - tau_1s(j+1)) / dSOC) * (x(1) - 0.5) +tau_1s(j);
        R_2 = ((R_2s(j) - R_2s(j+1)) / dSOC) * (x(1) - 0.5) + R_2s(j);
        tau_2 = ((tau_2s(j) - tau_2s(j+1)) / dSOC) * (x(1) - 0.5) +tau_2s(j);
    elseif (x(1) < 0.3) && (x(1) >= 0.1)
        j = 4;
        R_0 = ((R_0s(j) - R_0s(j+1)) / dSOC) * (x(1) - 0.3) + R_0s(j);
        R_1 = ((R_1s(j) - R_1s(j+1)) / dSOC) * (x(1) - 0.3) + R_1s(j);
        tau_1 = ((tau_1s(j) - tau_1s(j+1)) / dSOC) * (x(1) - 0.3) +tau_1s(j);
        R_2 = ((R_2s(j) - R_2s(j+1)) / dSOC) * (x(1) - 0.3) + R_2s(j);
        tau_2 = ((tau_2s(j) - tau_2s(j+1)) / dSOC) * (x(1) - 0.3) +tau_2s(j);
    else

        R_0 = R_0s(5);
        R_1 = R_1s(5);
        tau_1 = tau_1s(5);
        R_2 = R_2s(5);
        tau_2 = tau_2s(5);
   
    end 

    
    u = -data.I_mA(i)/1000; %Define positve current as Discharging
    A = diag([1;exp(-dt(i)/tau_1);exp(-dt(i)/tau_2)]);
    B = [-dt(i)/(Q_nom*3600);R_1*(1-exp(-dt(i)/tau_1));R_2*(1-exp(-dt(i)/tau_2))];
    x = A*x + B*u; %Prediction
    Sigma = A*Sigma*A.' + R; %Uncertianty
    C = [(3*OCV_SOC.p1*x(1)^2 + 2*OCV_SOC.p2*x(1) + OCV_SOC.p3), ... 
        -1, -1]; %Jacobian
    K = Sigma*C.'*inv(C*Sigma*C.'+ Q); %Kalman Gain 
    V_predict = (OCV_SOC.p1*x(1)^3 + OCV_SOC.p2*x(1)^2 ...
        + OCV_SOC.p3*x(1) + OCV_SOC.p4) - x(2) - x(3) - (u/1000)*R_0; 
    x = x + K*(data.Ecell_V(i) - V_predict); %Update
    SOCs(i) = x(1); 
    Sigma = (eye(3) - K*C)*Sigma; %Updated Uncertianty

end     
end 