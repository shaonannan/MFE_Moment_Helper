function [VEs,VIs] = VQs(fE,fI,etaE,etaI,SEE,SIE,SEI,SII,DEE,DIE,DEI,DII,LEE,LIE,mE,mI,NMDAE,NMDAI,vL,NE,NI,gL)
                  
%% Assume that vL = 0; gL = 0.005
% if (gL ==0)
%     warning('gL can not be 0 in VQs function');
% end
%% Self-connection
VEs(1,1) = (vL*gL + etaE(1,1)*fE(1,1) + NE*mE(1,1)*SEE - NI*mI(1,1)*SEI); 
VIs(1,1) = (vL*gL + etaI(1,1)*fI(1,1) + NE*mE(1,1)*SIE - NI*mI(1,1)*SII);
VEs(2,1) = (vL*gL + etaE(2,1)*fE(2,1) + NE*mE(2,1)*SEE - NI*mI(2,1)*SEI); 
VIs(2,1) = (vL*gL + etaI(2,1)*fI(2,1) + NE*mE(2,1)*SIE - NI*mI(2,1)*SII);
%% Cross-connection
VEs(1,1) = (VEs(1,1) + NE*mE(2)*DEE - NI*mI(2,1)*DEI); 
VIs(1,1) = (VIs(1,1) + NE*mE(2)*DIE - NI*mI(2,1)*DII);
VEs(2,1) = (VEs(2,1) + NE*mE(1)*DEE - NI*mI(1,1)*DEI); 
VIs(2,1) = (VIs(2,1) + NE*mE(1)*DIE - NI*mI(1,1)*DII);
%% Long-range NMDA type connection
VEs(1,1) = (VEs(1,1) + NMDAE(1,1)*LEE);
VIs(1,1) = (VIs(1,1) + NMDAI(1,1)*LIE);
VEs(2,1) = (VEs(2,1) + NMDAE(2,1)*LEE);
VIs(2,1) = (VIs(2,1) + NMDAI(2,1)*LIE);
%% TAU_M/GL
VEs = VEs/gL; VIs = VIs/gL;
