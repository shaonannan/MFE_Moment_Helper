function [NMDAE,NMDAI] = LRNMDA(NMDAE,NMDAI,mE,mI,gL,NE,NI,tau_r,tau_d)
dt = 0.1; trise = tau_r; tdamp = tau_d;
tr  = dt/trise;
etr = exp(-tr);
td  = dt/tdamp;
etd = exp(-td);
cst = 1.0/(tdamp-trise)*(etd-etr) ;   
%% Assume that vL = 0; gL = 0.005
HNMDAE = zeros(2,1); HNMDAI = zeros(2,1);
%% Hnmda
HNMDAE(1) = mE(2)*NE; HNMDAE(2) = mE(1)*NE;
HNMDAI(1) = mE(2)*NE; HNMDAI(2) = mE(1)*NE;
%% Inmda
NMDAE = NMDAE * etd + HNMDAE * cst;
NMDAI = NMDAI * etd + HNMDAI * cst;
