clear;clc
v_o = linspace(-1.0,1.0,200);

S_EE = 0.308*1.05/100;
S_IE = 0.308*1.05/100;
S_EI = 0.0363*1.05/100;
S_II = 0.0363*1.05/100;
D_EE = 0.03/100;
D_IE = 0.03/100;
D_EI = 0.150/100;
D_II = 0.150/100;
L_EE = 0%0.0102*1.28/100;
L_IE = 0%0.0612*1.28/100;
NE = 100;NI = 100;
LE = [4,0]; LI = [6,0];
load 20180903184415rhov.mat
rho_o = rhovE; rho_ln_o = rhovI;
NMDAE = zeros(2,1); NMDAI = zeros(2,1);
mEY = 3.11*0.76; mIY = 2.936*0.76;
DIY = 0.013145/0.76;  DEY = 0.013160/0.76;
mE = zeros(2,1); mI = zeros(2,1);
mE(1) = 0.01; mI(1) = 0.02;
for i = 1:1:100000
[mE_ra,mI_ra,rho_o,rho_ln_o,NMDAE,NMDAI] = Transient_2D_SDE_Solver_v1(v_o,rho_o,rho_ln_o,mE,mI,NMDAE,NMDAI,L_EE,L_IE,LE,LI,S_EE,S_IE,S_EI,S_II,D_EE,D_IE,D_EI,D_II,...
    mEY,mIY,DEY,DIY,NE,NI);
if mod(i,10) == 0
mE_ra,mI_ra,
end
mE = mE_ra; mI = mI_ra;
% pause;
end
