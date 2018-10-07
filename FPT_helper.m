%%% >>>>>>>>>>>>>>>>>>>>> Init >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
FlagInit = 1;format long;
if FlagInit >0
%% Hyper Parameters 
NE = 100;NI = 100; NPATCH = 2;
vL =0; gL=0.05;
Final_time = 1;dt = 0.1; TMAX = 1000;
SEE = 0.308*1.05;SEI = 0.0363*1.05;SIE = 0.308*1.05; SII = 0.0363*1.05; 
DEE = 0.03; DIE = 0.03; DEI = 0.150; DII = 0.150; 
LEE = 0.0102*1.28; LIE = 0.0612*1.28;
SEE = SEE/NE; SEI = SEI/NI; SIE = SIE/NE; SII = SII/NI;
DEE = DEE/NE; DIE = DIE/NE; DEI = DEI/NI; DII = DII/NI;
LEE = LEE/NE; LIE = LIE/NE;

etaE = 3.11*0.76 * ones(2,1); fE = 0.01316/0.76* ones(2,1);etaI = 2.36*0.76* ones(2,1);fI = 0.013145/0.76* ones(2,1); 
%% MFE and Vbins Vedges 
V_start = -1; V_end = 1; N_divide =201; V = linspace(V_start,V_end,N_divide); V = 0.5*V(2:end) + 0.5*V(1:end-1);h = V(2) - V(1);V = V';
step = 0; First_firing_num_neruons = 2; P1_used = 1;
counter = 1; t =0; tau_r = 2.0;tau_d = 128.0;
vT_idx = length(V); vT = 1.0; VR =0;
%% Load Data
load('20180903184415rhov.mat')
rE = rhovE; rI = rhovI;
%% To define some index, it is easy to find some index for V>=VT and and the index for v = VR in the coming code
rE  = rE./(h*sum(rE,2));rI = rI./(h*sum(rI,2));
vbarE = (V(2)-V(1))*(V'*rE');wbarE = (V(2)-V(1))*(V'.^2 *rE');vbar3E = (V(2)-V(1))*(V'.^3 *rE');vbar4E = (V(2)-V(1))*(V'.^4 *rE');
vbarI = (V(2)-V(1))*(V'*rI');wbarI = (V(2)-V(1))*(V'.^2 *rI');vbar3I = (V(2)-V(1))*(V'.^3 *rI');vbar4I = (V(2)-V(1))*(V'.^4 *rI');
%% transpose
vbarE = vbarE';vbarI = vbarI';wbarE = wbarE';wbarI = wbarI';vbar3E = vbar3E';vbar3I = vbar3I';vbar4E = vbar4E';vbar4I = vbar4I';
% % vbarE,vbarI,wbarE,wbarI,vbar3E,vbar3I,vbar4E,vbar4I,
% % pause;

%% Hold all; plot(V,source),hold on
NMDAE = zeros(2,1); NMDAI = zeros(2,1); mE = zeros(2,1); mI = zeros(2,1);
VEs = zeros(2,1); VIs = zeros(2,1); DE = zeros(2,1); DI = zeros(2,1);
PEq = zeros(2,length(V)); PIq = zeros(2,length(V)); sumE = zeros(2,1); sumI = zeros(2,1);
La0 = zeros(2,3); LaI0 = zeros(2,3); La1 = zeros(2,3); LaI1 = zeros(2,3);
% % % [VEs,VIs]  = VQs(fE,fI,etaE,etaI,SEE,SIE,SEI,SII,DEE,DIE,DEI,DII,LEE,LIE,mE,mI,NMDAE,NMDAI,vL,NE,NI,gL);
% % % [DE,DI]    = DEDI(fE,fI,etaE,etaI,SEE,SIE,SEI,SII,DEE,DIE,DEI,DII,mE,mI,NE,NI,gL);
% % % VEs,VIs,DE,DI,
% % % pause;
for idxPop = 1:1:2
% % % [PEq(idxPop),sumE(idxPop) ]= rho_EQ(VEs(idxPop),DE(idxPop),V);
% % % [PIq(idxPop),sumI(idxPop)] = rho_EQ(VIs(idxPop),DI(idxPop),V);
gammaE     = [1,vbarE(idxPop),wbarE(idxPop),vbar3E(idxPop),vbar4E(idxPop)];
gammaI     = [1,vbarI(idxPop),wbarI(idxPop),vbar3I(idxPop),vbar4I(idxPop)];
fiE = gammaE';     fiI = gammaI';
moment2 = 1;   options = optimset('TolFun',1e-9,'GradObj','on');
if moment2
F = fiE(2:3);FI = fiI(2:3);
else
F = fiE(2:end)*0;FI = fiI(2:end)*0;
end
if moment2 
La0(idxPop,:)  = gammaE(1:3)';
LaI0(idxPop,:) = gammaI(1:3)';
else
La0(idxPop,:)  = gammaE(1:end)';
LaI0(idxPop,:) = gammaI(1:end)';
end
end
nbins = 200; nbinp = nbins+1; dV = (vT-(VR-vT))/nbins; Vedges = linspace(VR-vT,vT,nbinp); Vedges = Vedges';Vbins = (Vedges(1:end-1)+Vedges(2:end))/2; 
j_source = ceil((nbins+1)/2); rho_source = zeros(nbins,1); rho_source(j_source)=1; 

N     = length(La0);
fin   = zeros(length(V),N);     fin(:,1)   = ones(size(V)); % fi0(x)=1
fin_s = zeros(length(Vbins),N); fin_s(:,1) = ones(size(Vbins)); % fi0(x)=1

 for n=2:N
      fin(:,n)   = V.*fin(:,n-1);
      fin_s(:,n) = Vbins.*fin_s(:,n-1);
 end
times = 10000; counter_firing_step = 0;
%%% >>>>>>>>>>>>> Recording Matrix >>>>>>>>>>>>>
plot_flag=0; 
dt_record_flag = 0; dtbin_record_flag = 1; 
tbinsize = 1.0; dtperbin = floor(tbinsize/dt); tbinsize=dtperbin*dt;
iteration_max = dtperbin * TMAX /tbinsize; 
%% *_ra for 
t_sum=0;    t_ra = zeros(iteration_max+1,1);
mE_ra = zeros(iteration_max+1,2);
mI_ra = zeros(iteration_max+1,2);
%% record each dt
if dt_record_flag;
rE_ra = zeros(2,nbins,iteration_max);
rI_ra = zeros(2,nbins,iteration_max);
else;%if dt_record_flag;
rE_ra = [];
rI_ra = [];
end;%if dt_record_flag;
%% record each time-bin
if dtbin_record_flag;
tbin_ra  = zeros(iteration_max/dtperbin,2);
mEbin_ra = zeros(iteration_max/dtperbin,2);
mIbin_ra = zeros(iteration_max/dtperbin,2);
xEbin_ra = zeros(iteration_max/dtperbin,2);
xIbin_ra = zeros(iteration_max/dtperbin,2);
rEbin_ra = zeros(2,nbins,iteration_max/dtperbin);
rIbin_ra = zeros(2,nbins,iteration_max/dtperbin);
VEavgbin_ra = zeros(iteration_max/dtperbin,2);
VEstdbin_ra = zeros(iteration_max/dtperbin,2);
VIavgbin_ra = zeros(iteration_max/dtperbin,2);
VIstdbin_ra = zeros(iteration_max/dtperbin,2);
P_MFEbin_ra = zeros(iteration_max/dtperbin,2);
else;%if dtbin_record_flag;
tbin_ra = [];
mEbin_ra = [];
mIbin_ra = [];
xEbin_ra = [];
xIbin_ra = [];
rEbin_ra = [];
rIbin_ra = [];
VEavgbin_ra = []; 
VEstdbin_ra = [];
VIavgbin_ra = [];
VIstdbin_ra = [];
P_MFEbin_ra = [];
end;%;%if dtbin_record_flag;
FlagInit = 0 ; %% Initiation finished
end

iteration = 1;
idx_kick = 2; idx_rec =1 ;N_neuron_fire_first=2;num_kick = 12;
while iteration<iteration_max;
t_ra(1+iteration) = t_sum;
t_sum
%% >>>>>>>>>> All about recording >>>>>>>>>> 
VIavg_ra(1+iteration,:) = vbarI;VIstd_ra(1+iteration,:) = sqrt(wbarI -vbarI.^2);
VEavg_ra(1+iteration,:) = vbarE;VEstd_ra(1+iteration,:) = sqrt(wbarE - vbarE.^2);
if dt_record_flag; rE_ra(:,:,1+iteration) = rE; rI_ra(:,:,1+iteration) = rI; end;
if dtbin_record_flag; 
dtbin_ij = 1 + floor(iteration/dtperbin);
tbin_ra(dtbin_ij) = tbin_ra(dtbin_ij) + dt*t_ra(1+iteration);
rEbin_ra(:,:,dtbin_ij)  = rEbin_ra(:,:,dtbin_ij) + dt*rE; 
rIbin_ra(:,:,dtbin_ij)  = rIbin_ra(:,:,dtbin_ij) + dt*rI; 
VEavgbin_ra(dtbin_ij,:) = VEavgbin_ra(dtbin_ij,:) + dt*VEavg_ra(1+iteration,:);
VEstdbin_ra(dtbin_ij,:) = VEstdbin_ra(dtbin_ij,:) + dt*VEstd_ra(1+iteration,:);
VIavgbin_ra(dtbin_ij,:) = VIavgbin_ra(dtbin_ij,:) + dt*VIavg_ra(1+iteration,:);
VIstdbin_ra(dtbin_ij,:) = VIstdbin_ra(dtbin_ij,:) + dt*VIstd_ra(1+iteration,:);
end;% if dtbin_record_flag;
%%% >>> >>> >>> >>> MFE >>> >>> >>> >>>
rho_o = rE;rho_ln_o = rI; 
for i = 1:1:2
rho_o(i,:) = rho_o(i,:)/sum(rho_o(i,:));
rho_ln_o(i,:) = rho_ln_o(i,:)/sum(rho_ln_o(i,:));
end

S_EE = SEE;S_IE = SIE;S_EI = SEI;S_II = SII;
v_o = linspace(-1.0,1.0,200);
%%% Calculate L_E,L_I
for i = 1:1:num_kick
[L_E,L_I,L_E_max,L_I_max,M1s_val,M2s_val,Ds_val,rE_o,rI_o,v_o] = Total_FPT_2D_SDE_Solver_old(v_o,rho_o(idx_kick,:),rho_ln_o(idx_kick,:),S_EE,S_IE,S_EI,S_II,NE,NI,N_neuron_fire_first);
%%% Refresh rhoV and NMDA
% S_EE,S_IE,S_EI,S_II,NE,NI,N_neuron_fire_first
LE = zeros(2,1);LI = zeros(2,1);
LE(idx_kick,1) = L_E+N_neuron_fire_first; LI(idx_kick,1) = L_I;
% % % LE,LI,
% % % pause;
for j = 1:1:2
mE_ra(iteration,j) = mE_ra(iteration,j) + LE(j,1)*dt;
mI_ra(iteration,j) = mI_ra(iteration,j) + LI(j,1)*dt;
end
[mEra,mIra,rho_o,rho_ln_o,NMDAE,NMDAI] = Transient_2D_SDE_Solver_v2(v_o,rho_o,rho_ln_o,NMDAE,NMDAI,LEE,LIE,LE,LI,S_EE,S_IE,S_EI,S_II,DEE,DIE,DEI,DII,...
    etaE(1),etaI(1),fE(1),fI(1),NE,NI);
figure(1);
subplot(2,1,1);
hold on;
plot(v_o,rho_o(1,:),'r',v_o,rho_ln_o(1,:),'b');
ylim([0,0.04]);
subplot(2,1,2);
hold on;
plot(v_o,rho_o(2,:),'r',v_o,rho_ln_o(2,:),'b');
ylim([0,0.04]);
pause(0.05);
end
rE = rho_o/(V(2)-V(1)); rI = rho_ln_o/(V(2)-V(1));
%%% >>> >>> >>> >>> END >>> >>> >>> >>>
%% Choose the population to fire(MFE)
% % % rho_o = rE;rho_ln_o = rI; 
% % % for i = 1:1:2
% % % rho_o(i,:) = rho_o(i,:)/sum(rho_o(i,:));
% % % rho_ln_o(i,:) = rho_ln_o(i,:)/sum(rho_ln_o(i,:));
% % % end
% % % 
% % % S_EE = SEE;S_IE = SIE;S_EI = SEI;S_II = SII;
% % % v_o = linspace(-1.0,1.0,200);
% % % %%% Calculate L_E,L_I
% % % for i = 1:1:num_kick
% % % [L_E,L_I,L_E_max,L_I_max,M1s_val,M2s_val,Ds_val,rE_o,rI_o,v_o] = Total_FPT_2D_SDE_Solver_old(v_o,rho_o(idx_kick,:),rho_ln_o(idx_kick,:),S_EE,S_IE,S_EI,S_II,NE,NI,N_neuron_fire_first);
% % % %%% Refresh rhoV and NMDA
% % % % S_EE,S_IE,S_EI,S_II,NE,NI,N_neuron_fire_first
% % % LE = zeros(2,1);LI = zeros(2,1);
% % % LE(idx_kick) = L_E+N_neuron_fire_first; LI(idx_kick) = L_I;
% % % for j = 1:1:2
% % % mE_ra(iteration,j) = mE_ra(iteration,j) + LE(j)*dt;mI_ra(iteration,j) = mI_ra(iteration,j) + LI(j)*dt;
% % % end
% % % [mEra,mIra,rho_o,rho_ln_o,NMDAE,NMDAI] = Transient_2D_SDE_Solver_v1(v_o,rho_o,rho_ln_o,NMDAE,NMDAI,LEE,LIE,LE,LI,S_EE,S_IE,S_EI,S_II,DEE,DIE,DEI,DII,...
% % %     etaE(1),etaI(1),fE(1),fI(1),NE,NI);
% % % figure(1);
% % % subplot(2,1,1);
% % % hold on;
% % % plot(v_o,rho_o(1,:),'r',v_o,rho_ln_o(1,:),'b');
% % % ylim([0,0.04]);
% % % subplot(2,1,2);
% % % hold on;
% % % plot(v_o,rho_o(2,:),'r',v_o,rho_ln_o(2,:),'b');
% % % ylim([0,0.04]);
% % % pause(0.1);
% % % end
% % % rE = rho_o/(V(2)-V(1)); rI = rho_ln_o/(V(2)-V(1));
%% Check Difference
[VEs,VIs]     = VQs(fE,fI,etaE,etaI,SEE,SIE,SEI,SII,DEE,DIE,DEI,DII,LEE,LIE,mE,mI,NMDAE,NMDAI,vL,NE,NI,gL);
DiffVEs = abs(VEs(2)-VEs(1));
DiffVEs = DiffVEs/max(VEs);
mE = zeros(2,1);mI = zeros(2,1);
flagDiff = 1; FirstofAll = 1;
while flagDiff >0
if FirstofAll >0
%% To define some index, it is easy to find some index for V>=VT and and the index for v = VR in the coming code
% rE  = rE./(h*sum(rE,2));rI = rI./(h*sum(rI,2));
vbarE = (V(2)-V(1))*(V'*rE');wbarE = (V(2)-V(1))*(V'.^2 *rE');vbar3E = (V(2)-V(1))*(V'.^3 *rE');vbar4E = (V(2)-V(1))*(V'.^4 *rE');
vbarI = (V(2)-V(1))*(V'*rI');wbarI = (V(2)-V(1))*(V'.^2 *rI');vbar3I = (V(2)-V(1))*(V'.^3 *rI');vbar4I = (V(2)-V(1))*(V'.^4 *rI');
%% transpose
vbarE = vbarE';vbarI = vbarI';wbarE = wbarE';wbarI = wbarI';vbar3E = vbar3E';vbar3I = vbar3I';vbar4E = vbar4E';vbar4I = vbar4I';
%% Hold all; plot(V,source),hold on
mE = zeros(2,1); mI = zeros(2,1);
PEq = zeros(2,length(V)); PIq = zeros(2,length(V)); sumE = zeros(2,1); sumI = zeros(2,1);
La0 = zeros(2,3); LaI0 = zeros(2,3); La1 = zeros(2,3); LaI1 = zeros(2,3);
for idxPop = 1:1:2
gammaE     = [1,vbarE(idxPop),wbarE(idxPop),vbar3E(idxPop),vbar4E(idxPop)];
gammaI     = [1,vbarI(idxPop),wbarI(idxPop),vbar3I(idxPop),vbar4I(idxPop)];
fiE = gammaE';     fiI = gammaI';
moment2 = 1;   options = optimset('TolFun',1e-9,'GradObj','on');
if moment2
F = fiE(2:3);FI = fiI(2:3);
else
F = fiE(2:end)*0;FI = fiI(2:end)*0;
end
if moment2 
La0(idxPop,:)  = gammaE(1:3)';
LaI0(idxPop,:) = gammaI(1:3)';
else
La0(idxPop,:)  = gammaE(1:end)';
LaI0(idxPop,:) = gammaI(1:end)';
end
end
FirstofAll = 0;
end

[NMDAE,NMDAI] = LRNMDA(NMDAE,NMDAI,mE,mI,gL,NE,NI,tau_r,tau_d);
[VEs,VIs]     = VQs(fE,fI,etaE,etaI,SEE,SIE,SEI,SII,DEE,DIE,DEI,DII,LEE,LIE,mE,mI,NMDAE,NMDAI,vL,NE,NI,gL);
[DE,DI]       = DEDI(fE,fI,etaE,etaI,SEE,SIE,SEI,SII,DEE,DIE,DEI,DII,mE,mI,NE,NI,gL);
[vbarE,vbarI,wbarE,wbarI,vbar3E,vbar3I,vbar4E,vbar4I] = ...
    solveVbarWbar4(dt,vbarE, vbarI, wbarE, wbarI,vbar3E, vbar3I,vbar4E, vbar4I,VEs,VIs,DE, DI,mE,mI,gL);    

%% separately processing
for idxPop = 1:1:NPATCH
fiE = [1,vbarE(idxPop),wbarE(idxPop),vbar3E(idxPop),vbar4E(idxPop)]';
fiI = [1,vbarI(idxPop),wbarI(idxPop),vbar3I(idxPop),vbar4I(idxPop)]';
[PEq(idxPop,:),sumE(idxPop) ]= rho_EQ(VEs(idxPop),DE(idxPop),V);
[PIq(idxPop,:),sumI(idxPop)] = rho_EQ(VIs(idxPop),DI(idxPop),V);
if moment2
F = fiE(2:3);FI = fiI(2:3);
else
F = fiE(2:end);FI = fiI(2:end);
end
[La1temp]  = fminsearch(@(La)optfun(F,V,La,squeeze(PEq(idxPop,:))',fin,1),squeeze(La0(idxPop,:))',options);
La1(idxPop,:) = La1temp;
[LaI1temp] = fminsearch(@(LaII)optfun(FI,V,LaII,squeeze(PIq(idxPop,:))',fin,1),squeeze(LaI0(idxPop,:))',options);
LaI1(idxPop,:) = LaI1temp;
La0(idxPop,:)  = real(La1(idxPop,:));   LaI0(idxPop,:) = real(LaI1(idxPop,:));
     
RvE(idxPop,:) = reshape(PEq(idxPop,:),nbins,1).*exp(fin(:,1:N)*reshape(La1(idxPop,:),N,1)); % Calculate p(x)
RvI(idxPop,:) = reshape(PIq(idxPop,:),nbins,1).*exp(fin(:,1:N)*reshape(LaI1(idxPop,:),N,1)); % Calculate p(x)

RvE(idxPop,:) = RvE(idxPop,:)/((V(2)-V(1))*sum(RvE(idxPop,:)));
RvI(idxPop,:) = RvI(idxPop,:)/((V(2)-V(1))*sum(RvI(idxPop,:)));

mE(idxPop) = gL*sqrt(DE(idxPop)) *exp(sum(La1(idxPop,:)))/sumE(idxPop)/2;
mI(idxPop) = gL*sqrt(DI(idxPop)) *exp(sum(LaI1(idxPop,:)))/sumI(idxPop)/2;
end

t,mE,mI,
if length(find(isnan(mE)))>0
    pause;
end
for i = 1:1:2
mE_ra(iteration,i) = mE_ra(iteration,i) + mE(i); mI_ra(iteration,i) = mI_ra(iteration,i) + mI(i);
end  
PLOT_or_NOT = 1;
if step >5;
    t, mE,mI
    step = 0;
    if PLOT_or_NOT
    figure(111)
    subplot(211)
    plot(V(1:5:end),RvE(1,1:5:end),'m',V(1:5:end),RvE(2,1:5:end),'r');
    subplot(212)
    plot(V(1:5:end),RvI(1,1:5:end),'g',V(1:5:end),RvI(2,1:5:end),'b');
    drawnow;
    end
end
t = dt*counter;
counter = counter + 1; step = step +1;     
% mI_ra(1+iteration,:) = mI; mE_ra(1+iteration,:) = mE;
iteration = iteration + 1;
%% Check Difference
rE = RvE; rI = RvI;
[VEs,VIs] = VQs(fE,fI,etaE,etaI,SEE,SIE,SEI,SII,DEE,DIE,DEI,DII,LEE,LIE,mE,mI,NMDAE,NMDAI,vL,NE,NI,gL);
DiffVEs = abs(mE(2)-mE(1));
DiffVEs = DiffVEs/max(mE);
DiffVIs = abs(mI(2)-mI(1));
DoffVIs = DiffVIs/max(mI);
if (DiffVEs < 1e-2)%&(DiffVIs < 1e-3)
    flagDiff = flagDiff + 1;
end
if (DiffVEs < 1e-2) & (flagDiff >10)
    flagDiff = 0;
    alt_kick = idx_rec;
    alt_rec  = idx_kick;
    idx_kick = alt_kick;
    idx_rec  = alt_rec;
end
if mod(iteration,10) == 0
figure(11);
subplot(2,1,1);
hold on;
plot(2:1:iteration,mE_ra(1:iteration-1,1),'r',2:1:iteration,mI_ra(1:iteration-1,1),'b');
xlim([0,iteration_max]);
subplot(2,1,2);
hold on;
plot(2:1:iteration,mE_ra(1:iteration-1,2),'r',2:1:iteration,mI_ra(1:iteration-1,2),'b');
xlim([0,iteration_max]);
pause(0.05);
end
end
end %%  end of if (MFE_COMPUT_USED) %%
   
% % % if dtbin_record_flag; 
% % % dtbin_ij = 1 + floor(iteration/dtperbin);
% % % mEbin_ra(dtbin_ij) = mEbin_ra(dtbin_ij) + (1-MFE_flag)*mE_ra(1+iteration)*NE*dt + MFE_flag*LE_ra(MFE_num);
% % % mIbin_ra(dtbin_ij) = mIbin_ra(dtbin_ij) + (1-MFE_flag)*mI_ra(1+iteration)*NI*dt + MFE_flag*LI_ra(MFE_num);
% % % xEbin_ra(dtbin_ij) = xEbin_ra(dtbin_ij) + (1-MFE_flag)*psample(mE_ra(1+iteration)*NE*dt) + MFE_flag*LE_ra(MFE_num);
% % % xIbin_ra(dtbin_ij) = xIbin_ra(dtbin_ij) + (1-MFE_flag)*psample(mI_ra(1+iteration)*NI*dt) + MFE_flag*LI_ra(MFE_num);
% % % P_MFEbin_ra(dtbin_ij) = P_MFEbin_ra(dtbin_ij) + P_MFE_ra(1+iteration)*dt;
% % % end;% if dtbin_record_flag;
% % % 
% % % if plot_flag;
% % % subplot(2,2,1); hold on; if mod(iteration,512)==0; plot(Vbins,rE,'r-',Vbins,rI,'b-'); end; hold off;
% % % end;%if plot_flag;
% % % 
% % % t_sum = t_sum+dt;
% % % counter_firing_step = 1+counter_firing_step;
% % % sum_mE = mE + sum_mE;sum_mI = mI + sum_mI;
% % %     
% % % if(counter_firing_step*dt >= 1.0)
% % % if density_plot_or_not
% % % filename_rec2 =  [catalog_head1,catalog_str,var_fn2,'_',num2str(density_record_counter),'.txt'];
% % % fid2 = fopen(filename_rec2, 'wt');
% % % lnv = length(V);
% % % for i = 1:5:lnv
% % %     fprintf(fid2,'%16.5e    %16.5e    %16.5e\n', V(i), RvE(i), RvI(i));
% % % end
% % % fclose(fid2);
% % % density_record_counter = density_record_counter + 1;
% % % end
% % % mE_pn1 = sum_mE/counter_firing_step;
% % % mI_ln1 = sum_mI/counter_firing_step;
% % % fprintf(fid,'%16.4e    %16.4e     %16.4e   %16.0e     %16.0e\n', t, mE_pn1,mI_ln1,0,0);
% % % fprintf(fid1,'%16.4e    %16.4e     %16.4e   %16.4e     %16.4e    %16.4e    %16.4e     %16.4e   %16.4e\n', t, vbarE,wbarE,vbar3E,vbar4E,vbarI,wbarI,vbar3I,vbar4I);
% % % sum_mE = 0; sum_mI = 0.0; counter_firing_step = 0;
% % % end 
% % %    
% % % if plot_flag;
% % % xlim([Vedges(1) Vedges(end)]);
% % % subplot(2,2,2); plot(Vbins,rE,'r.-',Vbins,rI,'b.-');
% % % subplot(2,2,[3 4]); plot(t_ra,mE_ra,'r.-',t_ra,mI_ra,'b.-',t_ra,P_MFE_ra,'g.-');
% % % %if plot_flag; print('-depsc',sprintf('%s_FIG_A.eps',filename_base)); end;
% % % end;%if plot_flag;
% % % 
% % % plot_flag=0;
% % % if plot_flag;
% % % figure;
% % % plot(Vbins,rE,'r-',Vbins,rI,'b-');
% % % xlim([Vedges(j_source) Vedges(end)]);
% % % %if plot_flag; print('-depsc',sprintf('%s_FIG_B.eps',filename_base)); end;
% % % end;%if plot_flag;
% % % 
% % % plot_flag=0;
% % % if plot_flag & dt_record_flag;
% % % [rtmp,ctmp] = size(rE_ra); j_source = (rtmp+1)/2;
% % % subplot(4,1,1);imagesc(rE_ra(end:-1:j_source,:)*sparse(1:ctmp,1:ctmp,1./sum(rE_ra,1)),[0 4/1024]);
% % % [rtmp,ctmp] = size(rI_ra); j_source = (rtmp+1)/2;
% % % subplot(4,1,2);imagesc(rI_ra(end:-1:j_source,:)*sparse(1:ctmp,1:ctmp,1./sum(rI_ra,1)),[0 4/1024]);
% % % subplot(4,1,3); plot(t_ra,P_MFE_ra,'g-'); xlim([min(t_ra) max(t_ra)]);
% % % subplot(4,1,4); 
% % % hold on; 
% % % plot(t_ra,mE_ra,'r-',t_ra,mI_ra,'b-');
% % % l = line([Lt_ra;Lt_ra],[zeros(size(LE_ra));LE_ra]/NE); set(l,'LineWidth',1,'Color',[1 0 0]);
% % % l = line([Lt_ra;Lt_ra],[zeros(size(LI_ra));LI_ra]/NI); set(l,'LineWidth',3,'Color',[0 0 1]);
% % % xlim([min(t_ra) max(t_ra)]);
% % % hold off;
% % % end;%if plot_flag & dt_record_flag;
% % % 
% % % plot_flag=1;
% % % nsec_to_plot=1;
% % % if plot_flag & dtbin_record_flag;
% % % figure;
% % % tij=max(1,(TMAX-1024*nsec_to_plot)/tbinsize):TMAX/tbinsize;
% % % [rtmp,ctmp] = size(rEbin_ra); j_source = (rtmp+1)/2;
% % % subplot(5,1,1);imagesc(rEbin_ra(end:-1:j_source,tij)*sparse(1:length(tij),1:length(tij),1./sum(rEbin_ra(end:-1:j_source,tij),1)),[0 4/1024]); ylabel('VE');
% % % set(gca,'YTick',[1 j_source]);set(gca,'YTickLabel',{'VT','VR'});
% % % [rtmp,ctmp] = size(rIbin_ra); j_source = (rtmp+1)/2;
% % % subplot(5,1,2);imagesc(rIbin_ra(end:-1:j_source,tij)*sparse(1:length(tij),1:length(tij),1./sum(rIbin_ra(end:-1:j_source,tij),1)),[0 4/1024]); ylabel('VI');
% % % set(gca,'YTick',[1 j_source]);set(gca,'YTickLabel',{'VT','VR'});
% % % subplot(5,1,3);
% % % hold on;
% % % l=stairs(tbin_ra(tij),VEavgbin_ra(tij),'r-'); set(l,'LineWidth',2);
% % % l=stairs(tbin_ra(tij),VEavgbin_ra(tij)+VEstdbin_ra(tij),'r:'); set(l,'LineWidth',1);
% % % l=stairs(tbin_ra(tij),VEavgbin_ra(tij)-VEstdbin_ra(tij),'r:'); set(l,'LineWidth',1);
% % % l=stairs(tbin_ra(tij),VIavgbin_ra(tij),'b-'); set(l,'LineWidth',2);
% % % l=stairs(tbin_ra(tij),VIavgbin_ra(tij)+VIstdbin_ra(tij),'b:'); set(l,'LineWidth',1);
% % % l=stairs(tbin_ra(tij),VIavgbin_ra(tij)-VIstdbin_ra(tij),'b:'); set(l,'LineWidth',1);
% % % xlim([min(tbin_ra(tij)) max(tbin_ra(tij))]);
% % % ylim([VR vT]); ylabel('Vavg');
% % % set(gca,'YTick',[0 1]);set(gca,'YTickLabel',{'VR','VT'});
% % % hold off;
% % % subplot(5,1,4); stairs(tbin_ra(tij),P_MFEbin_ra(tij),'g-'); xlim([min(tbin_ra(tij)) max(tbin_ra(tij))]); ylabel('P/dt');
% % % subplot(5,1,5); hold on; stairs(tbin_ra(tij),xEbin_ra(tij),'r-'); stairs(tbin_ra(tij),xIbin_ra(tij),'b-'); xlim([min(tbin_ra(tij)) max(tbin_ra(tij))]); ylabel('#spk'); hold off; 
% % % print('-depsc',sprintf('%s_JW_full_FIG_A.eps',filename_base));
% % % print('-djpeg',sprintf('%s_JW_full_FIG_A.jpg',filename_base));
% % % end;%if plot_flag & dtbin_record_flag;
% % % 
% % % plot_flag=1;
% % % if plot_flag & dtbin_record_flag;
% % % figure;
% % % subplot(2,2,1); 
% % % hold on;
% % % rEtmp = mean(rEbin_ra,2); l=stairs(Vbins,rEtmp,'r-'); set(l,'LineWidth',2);
% % % rItmp = mean(rIbin_ra,2); l=stairs(Vbins,rItmp,'b-'); set(l,'LineWidth',2);
% % % xlim([Vbins(j_source) Vbins(end)]); %ylim([0 0.005]);
% % % hold off;
% % % subplot(2,2,2);
% % % hold on;
% % % hEtmp = hist(mEbin_ra,[1:1.5*NE]); l=plot(log2([1:1.5*NE]),log2(1+hEtmp),'r-'); set(l,'LineWidth',2);
% % % hItmp = hist(mIbin_ra,[1:1.5*NI]); l=plot(log2([1:1.5*NI]),log2(1+hItmp),'b-'); set(l,'LineWidth',2);
% % % %hEtmp = hist(mEbin_ra,[1:1.5*NE]); l=loglog(([1:1.5*NE]),(1+hEtmp),'r-'); set(l,'LineWidth',2);
% % % %hItmp = hist(mIbin_ra,[1:1.5*NI]); l=loglog(([1:1.5*NI]),(1+hItmp),'b-'); set(l,'LineWidth',2);
% % % hold off;
% % % subplot(2,2,[3,4]);
% % % h = spectrum.mtm; Fs = 1000*1/tbinsize; hpsd = psd(h,(VEavgbin_ra*NE + VIavgbin_ra*NI)/(NE+NI),'Fs',Fs); plot(hpsd);
% % % print('-depsc',sprintf('%s_JW_full_FIG_B.eps',filename_base));
% % % print('-djpeg',sprintf('%s_JW_full_FIG_B.jpg',filename_base));
% % % end;%if plot_flag & dtbin_record_flag;