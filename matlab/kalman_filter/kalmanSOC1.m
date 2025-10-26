function SOCs = kalmanSOC1(OCV_SOC,R_0s,R_1s,tau_1s,R_2s,tau_2s,Q_nom,mVoltage,mCurrent,t,paramsCov,x)
      
assert(isa(OCV_SOC, 'cfit'), 'Input must be 3rd Order SOC/OCV Curve Approximation');

%%Extended Kalman Filter

%x_predict = zeros(3,size(data,1));
%x_update = zeros(3,size(data,1));
%SOC_0 = 1;
Sigma = eye(3); %Initial State Uncertainty
SOCs = zeros(size(t,1),1);
Q = 0.1; %Measurement Variance
dt = zeros(size(t,1),1);

for i=1:size(t,1)-1
    dt(i+1) = t(i+1) - t(i,1); 
end 

for i=1:size(t,1)
    %Interpolate within SOC Ranges to Find Parameters 
    dSOC = 0.2;
    if  x(1) >= 0.9
        R_0 = R_0s(1);
        R_1 = R_1s(1);
        tau_1 = tau_1s(1);
        R_2 = R_2s(1);
        tau_2 = tau_2s(1);
        C = paramsCov(:,:,1);
    elseif (x(1) < 0.9) && (x(1) >= 0.7)
        j = 1;
        R_0 = ((R_0s(j) - R_0s(j+1)) / dSOC) * (x(1) - 0.9) + R_0s(j);
        R_1 = ((R_1s(j) - R_1s(j+1)) / dSOC) * (x(1) - 0.9) + R_1s(j);
        tau_1 = ((tau_1s(j) - tau_1s(j+1)) / dSOC) * (x(1) - 0.9) +tau_1s(j);
        R_2 = ((R_2s(j) - R_2s(j+1)) / dSOC) * (x(1) - 0.9) + R_2s(j);
        tau_2 = ((tau_2s(j) - tau_2s(j+1)) / dSOC) * (x(1) - 0.9) +tau_2s(j);
        C = paramsCov(:,:,2);
    elseif (x(1) < 0.7) && (x(1) >= 0.5)
        j = 2;
        R_0 = ((R_0s(j) - R_0s(j+1)) / dSOC) * (x(1) - 0.7) + R_0s(j);
        R_1 = ((R_1s(j) - R_1s(j+1)) / dSOC) * (x(1) - 0.7) + R_1s(j);
        tau_1 = ((tau_1s(j) - tau_1s(j+1)) / dSOC) * (x(1) - 0.7) +tau_1s(j);
        R_2 = ((R_2s(j) - R_2s(j+1)) / dSOC) * (x(1) - 0.7) + R_2s(j);
        tau_2 = ((tau_2s(j) - tau_2s(j+1)) / dSOC) * (x(1) - 0.7) +tau_2s(j);
        C = paramsCov(:,:,3);
    elseif (x(1) < 0.5) && (x(1) >= 0.3)
        j = 3;
        R_0 = ((R_0s(j) - R_0s(j+1)) / dSOC) * (x(1) - 0.5) + R_0s(j);
        R_1 = ((R_1s(j) - R_1s(j+1)) / dSOC) * (x(1) - 0.5) + R_1s(j);
        tau_1 = ((tau_1s(j) - tau_1s(j+1)) / dSOC) * (x(1) - 0.5) +tau_1s(j);
        R_2 = ((R_2s(j) - R_2s(j+1)) / dSOC) * (x(1) - 0.5) + R_2s(j);
        tau_2 = ((tau_2s(j) - tau_2s(j+1)) / dSOC) * (x(1) - 0.5) +tau_2s(j);
        C = paramsCov(:,:,4);
    elseif (x(1) < 0.3) && (x(1) >= 0.1)
        j = 4;
        R_0 = ((R_0s(j) - R_0s(j+1)) / dSOC) * (x(1) - 0.3) + R_0s(j);
        R_1 = ((R_1s(j) - R_1s(j+1)) / dSOC) * (x(1) - 0.3) + R_1s(j);
        tau_1 = ((tau_1s(j) - tau_1s(j+1)) / dSOC) * (x(1) - 0.3) +tau_1s(j);
        R_2 = ((R_2s(j) - R_2s(j+1)) / dSOC) * (x(1) - 0.3) + R_2s(j);
        tau_2 = ((tau_2s(j) - tau_2s(j+1)) / dSOC) * (x(1) - 0.3) +tau_2s(j);
        C = paramsCov(:,:,4);
    else
        R_0 = R_0s(5);
        R_1 = R_1s(5);
        tau_1 = tau_1s(5);
        R_2 = R_2s(5);
        tau_2 = tau_2s(5);
        C = paramsCov(:,:,5);
    end 

    u = -mCurrent(i); %Define positve current as Discharging (defined opposite in data)
    A = diag([1;exp(-dt(i)/tau_1);exp(-dt(i)/tau_2)]);
    B = [-dt(i)/(Q_nom*3600);R_1*(1-exp(-dt(i)/tau_1));R_2*(1-exp(-dt(i)/tau_2))];
    x = A*x + B*u; %Prediction
    H = jacobian_params(R_1,tau_1/R_1,R_2,tau_2/R_2,x(2),x(3),u,dt(i));
    R = H*C*H.';
    Sigma = A*Sigma*A.' + R; %Uncertianty
    C = [(3*OCV_SOC.p1*x(1)^2 + 2*OCV_SOC.p2*x(1) + OCV_SOC.p3), ... 
        -1, -1]; %Jacobian
    K = Sigma*C.'*inv(C*Sigma*C.'+ Q); %Kalman Gain 
    V_predict = (OCV_SOC.p1*x(1)^3 + OCV_SOC.p2*x(1)^2 ...
        + OCV_SOC.p3*x(1) + OCV_SOC.p4) - x(2) - x(3) - u*R_0; 
    x = x + K*(mVoltage(i) - V_predict); %Update
    SOCs(i) = x(1); 
    Sigma = (eye(3) - K*C)*Sigma; %Updated Uncertianty
      
end 