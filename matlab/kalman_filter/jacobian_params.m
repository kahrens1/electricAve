function H = jacobian_params(R_1,C_1,R_2,C_2,V1,V2,I,dt)

 dV1_dR1 = (((-dt*exp(-dt/R_1*C_1) / (R_1^2*C_1)))*V1) + I*(1 - exp(-dt/R_1*C_1) - ((dt/R_1*C_1)*exp(-dt/R_1*C_1)));
 dV1_dC1 = (((-dt*exp(-dt/R_1*C_1) / (R_1*C_1^2)))*V1) + I*(1 - exp(-dt/R_1*C_1) - ((dt/R_1*C_1)*exp(-dt/R_1*C_1)));
 dV2_dR2 = (((-dt*exp(-dt/R_2*C_2) / (R_2^2*C_2)))*V2) + I*(1 - exp(-dt/R_2*C_2) - ((dt/R_2*C_2)*exp(-dt/R_2*C_2)));
 dV2_dC2 = (((-dt*exp(-dt/R_2*C_2) / (R_2*C_2^2)))*V2) + I*(1 - exp(-dt/R_2*C_2) - ((dt/R_2*C_2)*exp(-dt/R_2*C_2)));

H = [0,  0,    0,   0,   0;
     0,dV1_dR1,dV1_dC1,0,0;
     0,0,0 dV2_dR2 dV2_dC2];

end 