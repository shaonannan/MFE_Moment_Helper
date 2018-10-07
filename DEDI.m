function [DE,DI] = DEDI(fE,fI,etaE,etaI,SEE,SIE,SEI,SII,DEE,DIE,DEI,DII,mE,mI,NE,NI,gL)
        
%% Assume that vL = 0; gL = 0.005
% if (NI ==0 | NE == 0)
%     warning('NI or NE can not be 0 in VQs function');
% end
%% Self-connection
DE(1,1) = (etaE(1,1)*fE(1,1).^2 + NE*mE(1,1)*SEE.^2 + NI*mI(1,1)*SEI.^2); 
DI(1,1) = (etaI(1,1)*fI(1,1).^2 + NE*mE(1,1)*SIE.^2 + NI*mI(1,1)*SII.^2);
DE(2,1) = (etaE(2,1)*fE(2,1).^2 + NE*mE(2,1)*SEE.^2 + NI*mI(2,1)*SEI.^2); 
DI(2,1) = (etaI(2,1)*fI(2,1).^2 + NE*mE(2,1)*SIE.^2 + NI*mI(2,1)*SII.^2);
%% Cross-connection
DE(1,1) = (DE(1,1) + NE*mE(2,1)*DEE.^2 + NI*mI(2,1)*DEI.^2); 
DI(1,1) = (DI(1,1) + NE*mE(2,1)*DIE.^2 + NI*mI(2,1)*DII.^2);
DE(2,1) = (DE(2,1) + NE*mE(1,1)*DEE.^2 + NI*mI(1,1)*DEI.^2); 
DI(2,1) = (DI(2,1) + NE*mE(1,1)*DIE.^2 + NI*mI(1,1)*DII.^2);
%% TAU_M/GL
DE = DE/gL; DI = DI/gL;