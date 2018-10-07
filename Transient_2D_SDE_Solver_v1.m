function [mE_ra,mI_ra,rho_o,rho_ln_o,NMDAE,NMDAI] = Transient_2D_SDE_Solver_v1(v_o,rho_o,rho_ln_o,NMDAE,NMDAI,L_EE,L_IE,LE,LI,S_EE,S_IE,S_EI,S_II,D_EE,D_IE,D_EI,D_II,...
    mEY,mIY,DEY,DIY,NE,NI)
%%% >>>>>> Parameters for Standard Master Equation Method >>>>>>>>>>
TAU_V=20.0;VT=1;VR=0;dt = 0.1;
% % % if isempty(dt); dt = 0.1; end;
% >>> Short-range Connections >>>
%% Self-connection
DII = S_II; DIE = S_IE; DEI = S_EI; DEE = S_EE;
%% Cross-connection
SEE = D_EE;SIE = D_IE; SEI = D_EI; SII = D_II;
% % % if isempty(DIY); DIY = 0.013145/0.76; end;
% % % if isempty(DEY); DEY = 0.013160/0.76; end;
%% External-connection
DIY = DIY; DEY = DEY;
nbins = length(v_o); nbinp = nbins + 1;
j_source = ceil((nbins+1)/2);
Vedges = linspace(VR-VT,VT,nbinp);
% >>> Preparing for EPSP/IPSP Kicks
L_flow   = get_L_flow(dt,Vedges); L_flow1= L_flow;
LEY_kick = get_L_kick(DEY,Vedges);  lEY_fire = 1-transpose(sum(LEY_kick,1)); lEY_undr = transpose(sum(LEY_kick,1));
LEY_kick = LEY_kick -eye(size(LEY_kick)); flux_to_zero_LEY = -sum(LEY_kick,1);
LEY_kick(j_source,:) = LEY_kick(j_source,:) + flux_to_zero_LEY;
% flux_to_zero_LEY,
% pause;
LEE_kick = get_L_kick(DEE,Vedges);  lEE_fire = 1-transpose(sum(LEE_kick,1)); lEE_undr = transpose(sum(LEE_kick,1));
LEE_kick = LEE_kick -eye(size(LEE_kick)); flux_to_zero_LEE = -sum(LEE_kick,1);
LEE_kick(j_source,:) = LEE_kick(j_source,:) + flux_to_zero_LEE;
% flux_to_zero_LEE
% pause;
LEI_kick = get_L_kick(-DEI,Vedges); lEI_fire = 1-transpose(sum(LEI_kick,1)); lEI_undr = transpose(sum(LEI_kick,1));
LEI_kick = LEI_kick -eye(size(LEI_kick)); missing = -sum(LEI_kick,1);
LEI_kick(1,:) = LEI_kick(1,:) + missing;flux_to_zero_LEI = zeros(size(missing));
% flux_to_zero_LEI,
% pause;
LIY_kick = get_L_kick(DIY,Vedges);  lIY_fire = 1-transpose(sum(LIY_kick,1)); lIY_undr = transpose(sum(LIY_kick,1));
LIY_kick = LIY_kick -eye(size(LIY_kick)); flux_to_zero_LIY = -sum(LIY_kick,1);
LIY_kick(j_source,:) = LIY_kick(j_source,:) + flux_to_zero_LIY;
% LIY_kick
% pause;
LIE_kick = get_L_kick(DIE,Vedges);  lIE_fire = 1-transpose(sum(LIE_kick,1)); lIE_undr = transpose(sum(LIE_kick,1));
LIE_kick = LIE_kick -eye(size(LIE_kick)); flux_to_zero_LIE = -sum(LIE_kick,1);
LIE_kick(j_source,:) = LIE_kick(j_source,:) + flux_to_zero_LIE;
% LIE_kick
% pause;
LII_kick = get_L_kick(-DII,Vedges); lII_fire = 1-transpose(sum(LII_kick,1)); lII_undr = transpose(sum(LII_kick,1));
LII_kick = LII_kick -eye(size(LII_kick)); missing = -sum(LII_kick,1);
LII_kick(1,:) = LII_kick(1,:) + missing;flux_to_zero_LII = zeros(size(missing));
% LII_kick
% pause;
HEE_kick = get_L_kick(L_EE,Vedges); HEE_fire = 1-transpose(sum(HEE_kick,1)); HEE_undr = transpose(sum(HEE_kick,1));
HEE_kick = HEE_kick -eye(size(HEE_kick)); flux_to_zero_HEE = -sum(HEE_kick,1);
HEE_kick(j_source,:) = HEE_kick(j_source,:) + flux_to_zero_HEE;
% HEE_kick
% pause;
HIE_kick = get_L_kick(L_IE,Vedges); HIE_fire = 1-transpose(sum(HIE_kick,1)); HIE_undr = transpose(sum(HIE_kick,1));
HIE_kick = HIE_kick -eye(size(HIE_kick)); flux_to_zero_HIE = -sum(HIE_kick,1);
HIE_kick(j_source,:) = HIE_kick(j_source,:) + flux_to_zero_HIE;
% HIE_kick
% pause;
SEE_kick = get_L_kick(SEE,Vedges);  sEE_fire = 1-transpose(sum(SEE_kick,1)); sEE_undr = transpose(sum(SEE_kick,1));
SEE_kick = SEE_kick -eye(size(SEE_kick)); flux_to_zero_SEE = -sum(SEE_kick,1);
SEE_kick(j_source,:) = SEE_kick(j_source,:) + flux_to_zero_SEE;
% SEE_kick
% pause;
SEI_kick = get_L_kick(-SEI,Vedges); sEI_fire = 1-transpose(sum(SEI_kick,1)); sEI_undr = transpose(sum(SEI_kick,1));
SEI_kick = SEI_kick -eye(size(SEI_kick)); missing = -sum(SEI_kick,1);
SEI_kick(1,:) = SEI_kick(1,:) + missing;flux_to_zero_SEI = zeros(size(missing));
% SEI_kick
% pause;
SIE_kick = get_L_kick(SIE,Vedges);  sIE_fire = 1-transpose(sum(SIE_kick,1)); sIE_undr = transpose(sum(SIE_kick,1));
SIE_kick = SIE_kick -eye(size(SIE_kick)); flux_to_zero_SIE = -sum(SIE_kick,1);
SIE_kick(j_source,:) = SIE_kick(j_source,:) + flux_to_zero_SIE;
% SIE_kick
% pause;
SII_kick = get_L_kick(-SII,Vedges); sII_fire = 1-transpose(sum(SII_kick,1)); sII_undr = transpose(sum(SII_kick,1));
SII_kick = SII_kick -eye(size(SII_kick)); missing = -sum(SII_kick,1);
SII_kick(1,:) = SII_kick(1,:) + missing;flux_to_zero_SII = zeros(size(missing));
% SII_kick
% pause;
% % % %% Self-connection
% % % %% LE[1] --> Pop[1], LE[2] --> Pop[2]
% % % mEI_net = zeros(2,1); mII_net = zeros(2,1);
% % % mIE_net = zeros(2,1); mEE_net = zeros(2,1);
% % % 
% % % mIE_net(1) = NE*mE(1); mIE_net(2) = NE*mE(2);
% % % mEE_net(1) = NE*mE(1); mEE_net(2) = NE*mE(2);
% % % mII_net(1) = NI*mI(1); mII_net(2) = NI*mI(2);
% % % mEI_net(1) = NI*mI(1); mEI_net(2) = NI*mI(2);
% % % %% Cross-connection
% % % %% LE[1] --> Pop[2], LE[2] --> Pop[1]
% % % cEI_net = zeros(2,1); cII_net = zeros(2,1);
% % % cIE_net = zeros(2,1); cEE_net = zeros(2,1);
% % % 
% % % cIE_net(1) = NE*mE(2); cIE_net(2) = NE*mE(1);
% % % cEE_net(1) = NE*mE(2); cEE_net(2) = NE*mE(1);
% % % cII_net(1) = NI*mI(2); cII_net(2) = NI*mI(1);
% % % cEI_net(1) = NI*mI(2); cEI_net(2) = NI*mI(1);%% Self-connection
%% LE[1] --> Pop[1], LE[2] --> Pop[2]
mEI_net = zeros(2,1); mII_net = zeros(2,1);
mIE_net = zeros(2,1); mEE_net = zeros(2,1);

mIE_net(1) = LE(1)/dt; mIE_net(2) = LE(2)/dt;
% mIE_net,
% pause;
mEE_net(1) = LE(1)/dt; mEE_net(2) = LE(2)/dt;
mII_net(1) = LI(1)/dt; mII_net(2) = LI(2)/dt;
mEI_net(1) = LI(1)/dt; mEI_net(2) = LI(2)/dt;
%% Cross-connection
%% LE[1] --> Pop[2], LE[2] --> Pop[1]
cEI_net = zeros(2,1); cII_net = zeros(2,1);
cIE_net = zeros(2,1); cEE_net = zeros(2,1);

cIE_net(1,1) = LE(2,1)/dt; cIE_net(2,1) = LE(1)/dt;
cEE_net(1,1) = LE(2,1)/dt; cEE_net(2,1) = LE(1)/dt;
cII_net(1,1) = LI(2,1)/dt; cII_net(2,1) = LI(1)/dt;
cEI_net(1,1) = LI(2,1)/dt; cEI_net(2,1) = LI(1)/dt;

mE_ra = zeros(2,1); mI_ra = zeros(2,1);
%%% >>> Parameters >>>
deltat = dt; trise  = 2.0; tdamp  = 128.0;
tr   = deltat/trise; etr  = exp(-tr);
td   = deltat/tdamp; etd  = exp(-td);
cst  = 1.0/(tdamp - trise)*(etd - etr) ;
%%% >>> Long-range Connections >>>
HNMDAE = zeros(2,1); HNMDAI = zeros(2,1);
for i = 1:1:2
HNMDAE(i,1) = cEE_net(i,1); % mE_ra;
NMDAE(i,1)  = NMDAE(i,1) * etd + HNMDAE(i,1) * cst;
HNMDAI(i,1) = cIE_net(i,1); % mE_ra;
NMDAI(i,1)  = NMDAI(i,1) * etd + HNMDAI(i,1) * cst;
end

% % % rE = squeeze(rho_o(i,:));rE = rE';
% % % rI = squeeze(rho_ln_o(i,:));rI = rI';
% % % mEE_LR    = dt*NMDAE(i) * dot(HEE_fire,rE);
% % % mIE_LR    = dt*NMDAI(i) * dot(HIE_fire,rI);

% >>> Iteration >>>
% >>> For Excitatory Subpopulation >>>
for i = 1:1:2
rE = squeeze(rho_o(i,:));rE = rE';
rE       = L_flow*rE;
Vec_fire = zeros(size(flux_to_zero_LEY));
JmatE    = zeros(size(LEY_kick));
Vec_fire = Vec_fire + mEY*flux_to_zero_LEY ;
Vec_fire = Vec_fire + mEE_net(i)*flux_to_zero_LEE;
Vec_fire = Vec_fire + mEI_net(i)*flux_to_zero_LEI;
JmatE    = JmatE + mEY * LEY_kick ;
JmatE    = JmatE + mEE_net(i) * LEE_kick ;
JmatE    = JmatE + mEI_net(i) * LEI_kick ;
% VzerosE  = 
% mE_ra,
% pause;
% rEY_fire,rEE_fire,rEI_fire,
% size(JmatE),size(rE),size(Vec_fire)
rE       = expm(JmatE*dt)*rE;
mE_ra    = Vec_fire*rE;
% mE_ra,
% pause;
rho_o(i,:) = rE;
end

% >>> For Inhibitory Subpopulation >>>
for i = 1:1:2
rI = squeeze(rho_ln_o(i,:));rI = rI';
rI       = L_flow*rI;
Vec_fire = zeros(size(flux_to_zero_LIY));
JmatI    = zeros(size(LIY_kick));
Vec_fire = Vec_fire + mIY*flux_to_zero_LIY ;
Vec_fire = Vec_fire + mIE_net(i)*flux_to_zero_LIE;
Vec_fire = Vec_fire + mII_net(i)*flux_to_zero_LII;
JmatI    = JmatI + mIY * LIY_kick ;
JmatI    = JmatI + mIE_net(i) * LIE_kick ;
JmatI    = JmatI + mII_net(i) * LII_kick ;
% VzerosE  = 
% mE_ra,
% pause;
% rEY_fire,rEE_fire,rEI_fire,
% size(JmatI),size(rI),size(Vec_fire)
rI       = expm(JmatI*dt)*rI;
mI_ra    = Vec_fire*rI;
rho_ln_o(i,:) = rI;
end

%% >>> Cross-modification >>>
% >>> For Excitatory Subpopulation >>>
for i = 1:1:2
rE = squeeze(rho_o(i,:));rE = rE';
rE       = L_flow*rE;
Vec_fire = zeros(size(flux_to_zero_LEY));
JmatE    = zeros(size(LEY_kick));
Vec_fire = Vec_fire + mEY*flux_to_zero_LEY ;
Vec_fire = Vec_fire + cEE_net(i)*flux_to_zero_SEE;
Vec_fire = Vec_fire + cEI_net(i)*flux_to_zero_SEI;
Vec_fire = Vec_fire + NMDAE(i) * flux_to_zero_HEE;
JmatE    = JmatE + mEY * LEY_kick ;
JmatE    = JmatE + cEE_net(i) * SEE_kick ;
JmatE    = JmatE + cEI_net(i) * SEI_kick ;
JmatE    = JmatE + NMDAE(i) * HEE_kick ;
% VzerosE  = 
% mE_ra,
% pause;
% rEY_fire,rEE_fire,rEI_fire,
% size(JmatE),size(rE),size(Vec_fire)
rE       = expm(JmatE*dt)*rE;
mE_ra    = Vec_fire*rE;
% mE_ra,
% pause;
rho_o(i,:) = rE;
end

% >>> For Inhibitory Subpopulation >>>
for i = 1:1:2
rI = squeeze(rho_ln_o(i,:));rI = rI';
rI       = L_flow*rI;
Vec_fire = zeros(size(flux_to_zero_LIY));
JmatI    = zeros(size(LIY_kick));
Vec_fire = Vec_fire + mIY*flux_to_zero_LIY ;
Vec_fire = Vec_fire + cIE_net(i)*flux_to_zero_SIE;
Vec_fire = Vec_fire + cII_net(i)*flux_to_zero_SII;
Vec_fire = Vec_fire + NMDAI(i) * flux_to_zero_HIE;
JmatI    = JmatI + mIY * LIY_kick ;
JmatI    = JmatI + cIE_net(i) * SIE_kick ;
JmatI    = JmatI + cII_net(i) * SII_kick ;
JmatI    = JmatI + NMDAI(i) * HIE_kick ;
% mE_ra,
% pause;
% rEY_fire,rEE_fire,rEI_fire,
% size(JmatI),size(rI),size(Vec_fire)
rI       = expm(JmatI*dt)*rI;
mI_ra    = Vec_fire*rI;
rho_ln_o(i,:) = rI;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L_flow = get_L_flow(dt,Vedges);
%get L_flow -- the streaming term ;
TAU_V=20;VT=1;VR=0;
nbins = length(Vedges)-1; nbinp = nbins+1;
dV = (VT-(VR-VT))/nbins;
edt = exp(-dt/TAU_V); egt = exp(dt/TAU_V); 
row_ra=[];col_ra=[];val_ra=[];
for (j=1:nbins);
% bin edge j is at VR-VT + dV*(j-1); bin center j is at VR-VT + dV*(j-0.5) ;
% the edge index of V is j = (V-(VR-VT))/dV + 1 ;
% the edge center of V is j = (V-(VR-VT))/dV + 0.5 ;
% bin j with center Vbin(j) has edges [Vedges(j),Vedges(j+1)] ;
% over a timestep of length dt the mass in [Vpre Vpos] flows to lie in bin j, where ;
Vpre = Vedges(j)*egt;
Vpos = Vedges(j+1)*egt;
% The interval [Vpre Vpos] correponds to edges ;
jpre = (Vpre - (VR-VT))/dV + 1;
jpos = (Vpos - (VR-VT))/dV + 1;
% The interval [Vpre Vpos] sits within edge indices ;
jmin = floor(jpre);
jmax = ceil(jpos);
% and bin indices ;
bmin = jmin;
bmax = jmax-1;
% Because egt>1, jpos-jpre > dV, and so jmax-jmin>1 and bmax>bmin. ;
bvec = transpose(bmin:bmax);
% The weights associated with jmin and jmax are ;
wmin = (jmin+1)-jpre;
wmax = jpos-(jmax-1);
wvec = [wmin;ones(length(bvec)-2,1);wmax];
rvec = j*ones(length(bvec),1);
vij = find(bvec>0 & bvec<=nbins);
row_ra = [row_ra;rvec(vij)];
col_ra = [col_ra;bvec(vij)];
val_ra = [val_ra;wvec(vij)];
%if j==(nbins+1)/2; disp(sprintf('dV %0.6f Vpre %0.6f Vpos %0.6f jpre %0.6f jpos %0.6f jmin %d jmax %d bmin %d bmax %d ',dV,Vpre,Vpos,jpre,jpos,jmin,jmax,bmin,bmax)); end;
end;%for (j=1:nbins);
L_flow = sparse(row_ra,col_ra,val_ra,nbins,nbins);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L_kick = get_L_kick(kick_val,Vedges);
% get L_kick -- the kick term ;
TAU_V=20;VT=1;VR=0;
nbins = length(Vedges)-1; nbinp = nbins+1;
dV = (VT-(VR-VT))/nbins;
row_ra=[];col_ra=[];val_ra=[];
for (j=1:nbins);
% bin edge j is at VR-VT + dV*(j-1); bin center j is at VR-VT + dV*(j-0.5) ;
% the edge index of V is j = (V-(VR-VT))/dV + 1 ;
% the edge center of V is j = (V-(VR-VT))/dV + 0.5 ;
% bin j with center Vbin(j) has edges [Vedges(j),Vedges(j+1)] ;
% over a timestep of length dt the mass in [Vpre Vpos] is kicked to lie in bin j, where ;
Vpre = Vedges(j) - kick_val;
Vpos = Vedges(j+1) - kick_val;
% The interval [Vpre Vpos] correponds to edges ;
jpre = (Vpre - (VR-VT))/dV + 1;
jpos = (Vpos - (VR-VT))/dV + 1;
% The interval [Vpre Vpos] sits within edge indices ;
jmin = floor(jpre);
jmax = ceil(jpos);
% and bin indices ;
bmin = jmin;
bmax = jmax-1;
% because kicks do not compress, it is possible that jmax-jmin==1 and bmax==bmin. ;
bvec = transpose(bmin:bmax);
if length(bvec)>1;
% The weights associated with jmin and jmax are ;
wmin = (jmin+1)-jpre;
wmax = jpos-(jmax-1);
wvec = [wmin;ones(length(bvec)-2,1);wmax];
rvec = j*ones(length(bvec),1);
elseif length(bvec)==1;
% The weights associated with jmin and jmax are ;
wvec = [1];
rvec = j*ones(length(bvec),1);
end;%if length(bvec)>1;
vij = find(bvec>0 & bvec<=nbins);
row_ra = [row_ra;rvec(vij)];
col_ra = [col_ra;bvec(vij)];
val_ra = [val_ra;wvec(vij)];
end;%for (j=1:nbins);
L_kick = sparse(row_ra,col_ra,val_ra,nbins,nbins);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%