%%% Command to run
clear;clc
v_o = linspace(-1.0,1.0,200);

SEE  = 0.38*1.05/100;
S_IE = 0.38*1.05/100;
SEI  = 0.0363*1.05/100;
S_II = 0.0363*1.05/100;

L_EE = 0.0102*1.28/100;
L_IE = 0.0612*1.28/100;
NE = 100;NI = 100;
N_neuron_fire_first = 2;
% RECORDING
L_E = zeros(length(SEE),length(SEI));L_I = L_E;
L_E_max = L_E; L_I_max = L_I;
M1s_val = L_E; M2s_val = L_E; Ds_val = L_E;
load 20180903184415rhov.mat
rho_o = rhovE(1,:);rho_ln_o = rhovI(1,:);
rEp   = rhovE(2,:);rIp      = rhovI(2,:);
figure(1); 
subplot(2,1,1);
plot(v_o,rho_o,'r',v_o,rho_ln_o,'B');
subplot(2,1,2);
plot(v_o,rEp,'r',v_o,rIp,'B');
mEY = 3.11*0.76; mIY = 2.36*0.76;
recIEE_LR = 0; recIIE_LR = 0;
v_org = v_o;
for loop = 1:1:9
    loop,
for i = 1:1:length(SEE)
    for j = 1:1:length(SEI)
    S_EE = SEE(i); S_EI = SEI(j);
    S_IE = 1.0*S_EE; S_II = S_EI;
    [L_E,L_I,L_E_max,L_I_max,M1s_val,M2s_val,Ds_val,rho_o,rho_ln_o,v_o] = Total_FPT_2D_SDE_Solver(v_o,rho_o,rho_ln_o,S_EE,S_IE,S_EI,S_II,NE,NI,N_neuron_fire_first);
    v_o = (v_o)'; rho_o = (rho_o)';rho_ln_o = (rho_ln_o)';L_E = L_E + N_neuron_fire_first;  rEp = rEp';rIp = rIp';
    L_E,L_I
    pause;
    [mE_ra,mI_ra,rE,rI,rEp,rIp,recIEE_LR,recIIE_LR,mEp_ra,mIp_ra] = Total_Normal_2D_SDE_Solver(v_o,rho_o,rho_ln_o,S_EE,S_IE,S_EI,S_II,NE,NI,mEY,mIY,L_E,L_I,L_EE,L_IE,rEp,rIp,recIEE_LR,recIIE_LR);
    v_o = (v_o)'; rho_o = (rE)';rho_ln_o = (rI)'; rEp = rEp';rIp = rIp';
    figure(1);   
    subplot(2,1,1);%hold on;
    plot(v_org,rho_o,'.-r',v_o,rho_ln_o,'.-B'); 
    ylim([0.0,0.04]);
    subplot(2,1,2);%hold on;
    plot(v_o,rEp,'.-r',v_o,rIp,'.-B');
    ylim([0.0,0.04]);
    pause;    
    end
end
end



