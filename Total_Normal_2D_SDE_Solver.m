function [mE_ra,mI_ra,rE,rI,rEp,rIp,recIEE_LR,recIIE_LR,mEp_ra,mIp_ra] = Total_Normal_2D_SDE_Solver(v_o,rho_o,rho_ln_o,S_EE,S_IE,S_EI,S_II,NE,NI,mEY,mIY,LE,LI,L_EE,L_IE,rEp,rIp,recIEE_LR,recIIE_LR)
%%% >>>>>> Parameters for Standard Master Equation Method >>>>>>>>>>
TAU_V=20;VT=1;VR=0;
% % % if isempty(dt); dt = 0.1; end;
% >>> Short-range Connections >>>
DII = S_II; DIE = S_IE; DEI = S_EI; DEE = S_EE;

SEE = 0.03/NE;SIE = 0.03/NE; SEI = 0.150/NI; SII = 0.150/NI;

% % % if isempty(DIY); DIY = 0.013145/0.76; end;
% % % if isempty(DEY); DEY = 0.013160/0.76; end;
dt = 0.1;DIY = 0.013145/0.76; DEY = 0.013160/0.76;
nbins = length(v_o); nbinp = nbins + 1;
j_source = ceil((nbins+1)/2);
Vedges = linspace(VR-VT,VT,nbinp);
% >>> Preparing for EPSP/IPSP Kicks
L_flow = get_L_flow(dt,Vedges); L_flow1= L_flow;
LEY_kick = get_L_kick(DEY,Vedges);  lEY_fire = 1-transpose(sum(LEY_kick,1)); lEY_undr = transpose(sum(LEY_kick,1));
% lEY_fire(end-1:end),
% pause;
LEE_kick = get_L_kick(DEE,Vedges);  lEE_fire = 1-transpose(sum(LEE_kick,1)); lEE_undr = transpose(sum(LEE_kick,1));
% LEE_kick-eye(size(LEE_kick)),lEE_fire(end),
% pause;
LEI_kick = get_L_kick(-DEI,Vedges); lEI_fire = 1-transpose(sum(LEI_kick,1)); lEI_undr = transpose(sum(LEI_kick,1));
% LEI_kick-eye(size(LEE_kick)),lEI_fire(end),
% pause;
LIY_kick = get_L_kick(DIY,Vedges);  lIY_fire = 1-transpose(sum(LIY_kick,1)); lIY_undr = transpose(sum(LIY_kick,1));
LIE_kick = get_L_kick(DIE,Vedges);  lIE_fire = 1-transpose(sum(LIE_kick,1)); lIE_undr = transpose(sum(LIE_kick,1));
% LIE_kick-eye(size(LEE_kick)),lIE_fire(end),
% pause;
LII_kick = get_L_kick(-DII,Vedges); lII_fire = 1-transpose(sum(LII_kick,1)); lII_undr = transpose(sum(LII_kick,1));
% LII_kick-eye(size(LEE_kick)),lII_fire(end),
% pause;

HEE_kick = get_L_kick(L_EE,Vedges); HEE_fire = 1-transpose(sum(HEE_kick,1)); HEE_undr = transpose(sum(HEE_kick,1));
HIE_kick = get_L_kick(L_IE,Vedges); HIE_fire = 1-transpose(sum(HIE_kick,1)); HIE_undr = transpose(sum(HIE_kick,1));


SEE_kick = get_L_kick(SEE,Vedges);  sEE_fire = 1-transpose(sum(SEE_kick,1)); sEE_undr = transpose(sum(SEE_kick,1));
SEI_kick = get_L_kick(-SEI,Vedges); sEI_fire = 1-transpose(sum(SEI_kick,1)); sEI_undr = transpose(sum(SEI_kick,1));
SIE_kick = get_L_kick(SIE,Vedges);  sIE_fire = 1-transpose(sum(SIE_kick,1)); sIE_undr = transpose(sum(SIE_kick,1));
SII_kick = get_L_kick(-SII,Vedges); sII_fire = 1-transpose(sum(SII_kick,1)); sII_undr = transpose(sum(SII_kick,1));

rE = rho_o; mEE_net = LE/dt; mEI_net= LI/dt;
rI = rho_ln_o; mIE_net = LE/dt; mII_net= LI/dt;

% >>> Iteration >>>
% >>> For Excitatory Subpopulation >>>
% % % size(L_flow),size(rE),
rE       = L_flow*rE;
rEY_fire = dt*mEY*dot(lEY_fire,rE); 
rEE_fire = dt*mEE_net*dot(lEE_fire,rE);
rEI_fire = dt*mEI_net*dot(lEI_fire,rE);
mE_ra    = (rEY_fire + rEE_fire + rEI_fire)/dt;
% rEY_fire,rEE_fire,rEI_fire,
rE       = (1-dt*(mEY + mEE_net + mEI_net))*rE + dt*(mEY*LEY_kick*rE + mEE_net*LEE_kick*rE + mEI_net*LEI_kick*rE);
rE(j_source) = rE(j_source) + rEY_fire + rEE_fire + rEI_fire;
% % % Jmat1 = LEY_kick - eye(size(LEY_kick));
% % % Jmat2 = LEE_kick - eye(size(LEE_kick));
% % % Jmat3 = LEI_kick - eye(size(LEI_kick));
% % % rE_ = (rEo+dt*(mEY*Jmat1 + mEE_net*Jmat2+mEI_net*Jmat3)*rEo);%exp(dt*(mEY*Jmat1 + mEE_net*Jmat2+mEI_net*Jmat3))*rEo;
% % % rE_(j_source) = rE_(j_source) + rEY_fire + rEE_fire + rEI_fire;
% % % sum(rE_)
% % % pause;

% >>> For Inhibitory Subpopulation >>>
rI       = L_flow*rI;
rIY_fire = dt*mIY*dot(lIY_fire,rI); 
rIE_fire = dt*mIE_net*dot(lIE_fire,rI);
rII_fire = dt*mII_net*dot(lII_fire,rI);
mI_ra    = (rIY_fire + rIE_fire + rII_fire)/dt;
% rIY_fire, rIE_fire,rII_fire,
rI       = (1-dt*(mIY + mIE_net + mII_net))*rI + dt*(mIY*LIY_kick*rI + mIE_net*LIE_kick*rI + mII_net*LII_kick*rI);
rI(j_source) = rI(j_source) + rIY_fire + rIE_fire + rII_fire;

% % % mEE_net = (NE-1)*(rEY_fire + rEE_fire + rEI_fire)/dt;
% % % mIE_net = (NE-0)*(rEY_fire + rEE_fire + rEI_fire)/dt;
% % % %%% Two items above, both represent the capability of excitatory
% % % %%% subpopulation to have influence on others(both excitatory[self] and
% % % %%% inhibitory subpopulations.
% % % mEI_net = (NI-0)*(rIY_fire + rIE_fire + rII_fire)/dt;
% % % mII_net = (NI-1)*(rIY_fire + rIE_fire + rII_fire)/dt;
% % % %%% Two items above, both represent the capability of inhibitory
% % % %%% subpopulation ti have influence on others(both excitatory and
% % % %%% inhibitory

%{
%%% >>>>>>>>>>>>>>>>>> First Method~~~ >>>>>>>>>>>>>>>>>>>>>\%%%
%%% >>> Parameters >>>
deltat = dt; trise  = 2.0; tdamp  = 128.0;
tr   = deltat/trise; etr  = exp(-tr);
td   = deltat/tdamp; etd  = exp(-td);
cst  = 1.0/(tdamp - trise)*(etd - etr) ;
%%% >>> Long-range Connections >>>
recHEE_LR = NE * mE_ra;
recIEE_LR = recIEE_LR * etd + recHEE_LR * cst;
recHIE_LR = NE * mE_ra;
recIIE_LR = recIIE_LR * etd + recHIE_LR * cst;
mEE_LR    = dt*recIEE_LR * dot(HEE_fire,rEp);
mIE_LR    = dt*recIIE_LR * dot(HIE_fire,rIp);

% % % rEp    = (1-dt*(mEY + recIEE_LR))*rEp + dt*(mEY*LEY_kick*rEp + recIEE_LR*HEE_kick*rEp );
% % % rIp    = (1-dt*(mIY + recIIE_LR))*rIp + dt*(mIY*LIY_kick*rIp + recIIE_LR*HIE_kick*rIp );
rEp       = L_flow*rEp;
rEYp_fire = dt*mEY*dot(lEY_fire,rEp); 
rEEp_fire = dt*mE_ra*NE*dot(sEE_fire,rEp);
rEIp_fire = dt*mI_ra*NI*dot(sEI_fire,rEp);
mEp_ra    = (rEYp_fire + rEEp_fire + rEIp_fire + mEE_LR)/dt;
rEp       = (1-dt*(mEY + mE_ra*NE + mI_ra*NI +  + recIEE_LR))*rEp + dt*(mEY*LEY_kick*rEp + mE_ra*NE*SEE_kick*rEp + mI_ra*NI*SEI_kick*rEp + recIEE_LR*HEE_kick*rEp);
rEp(j_source) = rEp(j_source) + mEp_ra*dt;

% >>> For Inhibitory Subpopulation >>>
rIp       = L_flow*rIp;
rIYp_fire = dt*mIY*dot(lIY_fire,rIp); 
rIEp_fire = dt*mE_ra*NE*dot(sIE_fire,rIp);
rIIp_fire = dt*mI_ra*NI*dot(sII_fire,rIp);
mIp_ra    = (rIYp_fire + rIEp_fire + rIIp_fire + mIE_LR)/dt;
rIp       = (1-dt*(mIY + mE_ra*NE + mI_ra*NI + recIIE_LR))*rEp + dt*(mIY*LIY_kick*rIp + mE_ra*NE*SIE_kick*rIp + mI_ra*NI*SII_kick*rIp + recIIE_LR*HIE_kick*rIp );
rIp(j_source) = rIp(j_source) + mIp_ra*dt;
mEp_ra,mIp_ra,
%}

%{%}
%%% >>>>>>>>>>>>>>>>>> Second Method~~~ >>>>>>>>>>>>>>>>>>>>>\%%%
%%% >>> Parameters >>>
deltat = dt; trise  = 2.0; tdamp  = 128.0;
tr   = deltat/trise; etr  = exp(-tr);
td   = deltat/tdamp; etd  = exp(-td);
cst  = 1.0/(tdamp - trise)*(etd - etr) ;
%%% >>> Long-range Connections >>>
recHEE_LR = NE * mEE_net/NE; % mE_ra;
recIEE_LR = recIEE_LR * etd + recHEE_LR * cst;
recHIE_LR = NE * mIE_net/NE; % mE_ra;
recIIE_LR = recIIE_LR * etd + recHIE_LR * cst;
mEE_LR    = dt*recIEE_LR * dot(HEE_fire,rEp);
mIE_LR    = dt*recIIE_LR * dot(HIE_fire,rIp);

% % % rEp    = (1-dt*(mEY + recIEE_LR))*rEp + dt*(mEY*LEY_kick*rEp + recIEE_LR*HEE_kick*rEp );
% % % rIp    = (1-dt*(mIY + recIIE_LR))*rIp + dt*(mIY*LIY_kick*rIp + recIIE_LR*HIE_kick*rIp );
rEp       = L_flow*rEp;
rEYp_fire = dt*mEY*dot(lEY_fire,rEp); 
rEEp_fire = dt*mEE_net*dot(sEE_fire,rEp); % mE_ra*NE*dot(sEE_fire,rEp);
rEIp_fire = dt*mEI_net*dot(sEI_fire,rEp); % mI_ra*NI*dot(sEI_fire,rEp);
mEp_ra    = (rEYp_fire + rEEp_fire + rEIp_fire + mEE_LR)/dt;
rEp       = (1-dt*(mEY + mEE_net + mEI_net + recIEE_LR))*rEp + dt*(mEY*LEY_kick*rEp + mEE_net*SEE_kick*rEp + mEI_net*SEI_kick*rEp + recIEE_LR*HEE_kick*rEp);
% rEp       = (1-dt*(mEY + mE_ra*NE + mI_ra*NI +  + recIEE_LR))*rEp + dt*(mEY*LEY_kick*rEp + mE_ra*NE*SEE_kick*rEp + mI_ra*NI*SEI_kick*rEp + recIEE_LR*HEE_kick*rEp);
rEp(j_source) = rEp(j_source) + rEYp_fire + rEEp_fire + mEE_LR;%mEp_ra*dt;

% >>> For Inhibitory Subpopulation >>>
rIp       = L_flow*rIp;
rIYp_fire = dt*mIY*dot(lIY_fire,rIp); 
rIEp_fire = dt*mIE_net*dot(sIE_fire,rIp); % mE_ra*NE*dot(sIE_fire,rIp);
rIIp_fire = dt*mII_net*dot(sII_fire,rIp); % mI_ra*NI*dot(sII_fire,rIp);
mIp_ra    = (rIYp_fire + rIEp_fire + rIIp_fire + mIE_LR)/dt;
rIp       = (1-dt*(mIY + mIE_net + mII_net + recIIE_LR))*rIp + dt*(mIY*LIY_kick*rIp + mIE_net*SIE_kick*rIp + mII_net*SII_kick*rIp + recIIE_LR*HIE_kick*rIp );
rIp(j_source) = rIp(j_source) + rIYp_fire + rIEp_fire + mIE_LR;%mIp_ra*dt;
mEp_ra,mIp_ra,
%}


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%