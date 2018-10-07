function [LE,LI,LEmax,LImax,NMDAE,NMDAI,rho_o,rho_ln_o,v_o] = FPT_2D_SDE_Solver(v_o,rho_o,rho_ln_o,S_EE,S_IE,S_EI,S_II,NE,NI,NMDAE,NMDAI,N_neuron_fire_first)
N_sample = 10000;
vT =1;vI = -2/3;
tStep = 0.01; tStart = 0; tEnd = 1.0; 
tt = tStart:tStep:tEnd; dt = tStep;
alpha = S_EI/S_EE *NI/NE; NESEE = NE*S_EE;  
% h = v_o(2) - v_o(1);
% [cdf_pn, cdf_ln] = cdf(rho_o, rho_ln_o,h);
% G = cdf_pn + alpha*(1-cdf_ln);
% ll = 1-(vT - v_o)/NESEE+2/NE;
% figure; plot(v_o,G);
% hold all; plot(v_o,ll);

[rho,rho_ln, f, f_ln,v] = transform_from_rho_and_rho_ln(v_o,rho_o,rho_ln_o,S_EE,S_IE,vT,vI);
r = interp1(v,rho,tt); r_ln = interp1(v,rho_ln,tt);
% Re-scale
r = r/sum(r); r_ln = r_ln/sum(r_ln);
F = interp1(v,f,tt); F_ln = interp1(v,f_ln,tt);
% Re-scale
F = F/F(end); F_ln = F_ln/F_ln(end);
rdt = sqrt(dt);

count1 = 1; All_phi = [];
for jj = 1:N_sample; rNE = sqrt(NE*1); rNI = sqrt(NI*1);
    % flag, represents whether find the correct crossing value or not
    flag = 0; % Haven't find cross;
    counter = 1;t = counter*dt; phi10 = 0;phi20=0;
    
%     jj 
while (t<tEnd)
    rho_temp = r(counter);
    F_temp = F(counter); dw = randn(1,1);
    phi1 = evolution(phi10,F_temp,dt,rNE,rdt,dw,rho_temp);
    rho_ln_temp = r_ln(counter);
    F_ln_temp = F_ln(counter); dw = randn(1,1);
    phi2 = evolution(phi20,F_ln_temp,dt,rNI,rdt,dw,rho_ln_temp);
    t = counter*dt; 
    counter = counter+1; phi10 = phi1; phi20 = phi2;
    val = get_constrainted_plane_vale(phi1,phi2,t,alpha,F_temp,F_ln_temp,NESEE,NE,N_neuron_fire_first);
    
    if val>0
      All_phi(count1,1) = t; All_phi(count1,2) = phi1;All_phi(count1,3) = phi2;
%     All_phi(count1,1) = phi1+F_temp + alpha*(-F_ln_temp-phi2);
%     All_phi(count1,1) = F_temp + alpha*(-F_ln_temp);
%     All_phi(count1,1) = val;
%     All_phi(count1,1) = phi1;
    flag = 1;
    count1 = count1 + 1;
    break;
    end
end
end

M1s_val = 0;
M2s_val = 0;
Ds_val  = 0;

% figure;plot(tt(1:end-1),All_phi);
% ll = tt/NESEE -2/NE;
% hold on;
% plot(tt,ll)
% size(All_phi)
if count1>2
y = All_phi(:,1); x = All_phi(:,2);
% For mean-field theory
xdivided = 20/0.01; ydivided = 100; xs = -10; xe = 10.0; ys = 0; ye = 1;
mhist = hist2d([x'; y'],xdivided,ydivided,[xs,xe],[ys,ye]);
dx    = (xe - xs)/xdivided; dy = (ye-ys)/ydivided;
ycord = ys:dy:ye;
s_val = zeros(1,length(mhist(:,1)));
for jj = 1:length(mhist(:,1))
    s_val(jj) = dx * sum(mhist(jj,:));
end
sum_val = sum(s_val)*dy;
s_val = s_val/sum_val;
sval_max = max(s_val);
tval = tt(find(s_val == sval_max));
tval = tval(1);
% s_val,tval
ttIdx = find(tt<tval);
h = tt(2) -tt(1);
LE = sum(r(ttIdx))*NE;
LI = sum(r_ln(ttIdx))*NI;
% % % LE,
% % % LI,
% Calculate first-order Moment and second-order Moment
ySpace = linspace(ys,ye,ydivided);
M1s_val = sum(ySpace.*s_val)/length(s_val);
M2s_val = sum((ySpace.^2).*s_val)/length(s_val);
Ds_val  = M2s_val - M1s_val.^2;

% For max(fluctuation)
maxV_cross = min(y);
% % % maxV_cross
ttIdx = find(tt<maxV_cross);
h = tt(2) -tt(1);
LEmax = sum(r(ttIdx))*NE;
LImax = sum(r_ln(ttIdx))*NI;
% % % LEmax,
% % % LImax,
elseif count1>1
y = All_phi(:,1); x = All_phi(:,2);
ttIdx = find(tt<y);
h = tt(2) -tt(1);
LE = sum(r(ttIdx))*NE;
LI = sum(r_ln(ttIdx))*NI;
% % % LE,
% % % LI,

LEmax = LE;
LImax = LI;
% % % LEmax,
% % % LImax,
else
if flag <1
y = tt(end);
ttIdx = find(tt<y);
h = tt(2) -tt(1);
LE = sum(r(ttIdx))*NE;
LI = sum(r_ln(ttIdx))*NI;
% % % LE,
% % % LI,
    
LEmax = LE;
LImax = LI;
% % % LEmax,
% % % LImax,
end
end


% >>> Inverse, Preparing for Latter Processing >>>
[rho_o,rho_ln_o,v_o] = transform_to_rho_and_rho_ln(v,rho,rho_ln,S_EE,S_IE,vT,vI);
end %% end of sampled_exact_solution_used

function val=get_constrainted_plane_vale(phi1,phi2,t,alpha,F_temp,F_ln_temp,NESEE,NE,N_neuron_fire_first);
val = - phi1 + alpha*phi2 - F_temp + alpha*F_ln_temp + t/NESEE - N_neuron_fire_first/NE;
end

function phi = evolution(phi0,f,dt,r_N,rdt,dw,rho)
phi = phi0  - dt/r_N*((phi0.*rho)./(1.0 - f)) + rdt/r_N*(sqrt(rho).*dw);
% phi = phi0  - dt*((phi0.*rho)./(1.0 - f)) + 1/r_N*(sqrt(rho).*dw);
end

function [cdf_pn_out, cdf_ln_out] = cdf(rho, rho_ln,h)
for i = 1:length(rho)
    cdf_pn_out(i) = h*sum(rho(1:i));
    cdf_ln_out(i) = h*sum(rho_ln(1:i));
end
cdf_pn_out = cdf_pn_out/cdf_pn_out(end); %% this is the normalization to make sure than the probability is not great 1:
cdf_ln_out = cdf_ln_out/cdf_ln_out(end);
end



function [rho_out,rho_ln_out, f_out, f_ln_out,v] = transform_from_rho_and_rho_ln(v_o,rho_o,rho_ln_o,S_EE,S_IE,vT,vI);
% figure;
% subplot(4,1,1);
%  plot(v_o,rho_o,'r',v_o,rho_ln_o);
jj = 1:length(v_o);
v(jj) = v_o(end+1-jj); %% F: (vI,vT)--> (vT,vI)
v = vT - v; %% F:(vT,vI) --> (0,vT-vI)
rho_out(jj) = rho_o(end+1-jj); rho_ln(jj) = rho_ln_o(end+1-jj); %% transform the inverse direction for v in (0,vT-vI)
% subplot(4,1,2);
% plot(v,rho_out,'r',v,rho_ln);
%%% rho_ln_out(v) = SIE/SEE*rho_ln(SIE/SEE*v);

sIEE = S_IE/S_EE; s = sIEE*v; 
m_val = vT-vI;
for id = 1:length(s);
if(s(id)<=m_val)
rho_ln_out(id) = interp1(v,rho_ln,s(id));    
else
rho_ln_out(id) = 0;
end
end
rho_ln_out = rho_ln_out*sIEE;
% subplot(4,1,3);
% plot(v,rho_out,'r',v,rho_ln_out);
h = v(2) -v(1);
[f_out, f_ln_out] = cdf(rho_out, rho_ln_out,h);
% g(jj) = f_out(jj) + 1*(1-f_ln_out(jj));
% figure;plot(v,f_out,'r',v,f_ln_out,'g',v,g,'o')
% subplot(4,1,4);
% plot(v,f_out,'r',v,f_ln_out)
end


function [rho_out,rho_ln,v] = transform_to_rho_and_rho_ln(v_o,rho_o,rho_ln_o,S_EE,S_IE,vT,vI);
% inverse re-scale
sIEE = S_EE/S_IE; s = sIEE*v_o; 
rho_ln_out = zeros(size(rho_ln_o));
for id = 1:length(s);
rho_ln_out(id) = interp1(v_o,rho_ln_o,s(id));   
end
rho_ln_out = rho_ln_out*sIEE;

jj = 1:length(v_o);
v(jj) = v_o(end+1-jj); 
v = vT - v; 
rho_out(jj) = rho_o(end+1-jj); rho_ln(jj) = rho_ln_out(end+1-jj); 

end

function [t, rho,rho_ln] = get_initialized_rho_tt();
v_load = load('record1.txt');
t = v_load(:,1);rho = v_load(:,2);rho_ln = v_load(:,3);
end
