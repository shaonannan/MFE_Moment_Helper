function [vbarE1,vbarI1,wbarE1,wbarI1,vbar3E1,vbar3I1,vbar4E1,vbar4I1] ...
        = solveVbarWbar4(dt,vbarE, vbarI, wbarE, wbarI,vbar3E, vbar3I,vbar4E, vbar4I,VEs,VIs,DE, DI,mE,mI,gL)
%% Here assume thatn vT =1; so the coefficient of mE and mI is 1, we do not multiply them ...
dtgL =dt*gL;
vbarE1 = vbarE + dtgL*(- mE/gL - (vbarE - VEs));
vbarI1 = vbarI + dtgL*(- mI/gL - (vbarI - VIs));

wbarE1 = wbarE + dtgL*(- mE/gL - 2.0*(wbarE - VEs.*vbarE - 0.5*DE));
wbarI1 = wbarI + dtgL*(- mI/gL - 2.0*(wbarI - VIs.*vbarI - 0.5*DI) );

vbar3E1 = vbar3E + dtgL*(- mE/gL - 3.0*(vbar3E - VEs.*wbarE - DE.*vbarE));
vbar3I1 = vbar3I + dtgL*(- mI/gL - 3.0*(vbar3I - VIs.*wbarI - DI.*vbarI));

vbar4E1 = vbar4E + dtgL*(- mE/gL - 4.0*(vbar4E - VEs.*vbar3E - 1.5*DE.*wbarE));
vbar4I1 = vbar4I + dtgL*(- mI/gL - 4.0*(vbar4I - VIs.*vbar3I - 1.5*DI.*wbarI));