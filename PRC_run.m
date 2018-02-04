% script for computing the limit cycle trajectory and infinitesimal phase 
% response curve (iPRC) for an adaptive exponential integrate-and-fire 
% neuron, based on the method described in Ladenbauer et al., PLOS
% Computational Biology 2012; written by Josef Ladenbauer

tic 

%%% set parameters
C = 100;   % membrane capacitance [pF]
g_l = 10;   % leak conductance [nS]
E_l = -70;   % leak reversal potential [mV]
V_t = -50;   % spike threshold [mV]
delta_T = 2;   % slope factor [mV]
tau_w = 100;   % adaptation time constant [ms]
a = 15;   % subthreshold adaptation [nS]
b = 50;  % spike-triggered adaptation [pA]
V_r = -60;   % reset membrane voltage [mV] 
I_inj = 500;  % input current [pA]
pert_V = 0.1;  % voltage perturbation strength [mV]

w_init = 50;
T_max = 1000;  % max. period [ms]

%%% compute 
param = [C g_l E_l V_t delta_T tau_w a b V_r];
[t, lc] = aEIF_lcf(param,I_inj,V_r,w_init,T_max);  % limit cycle
if lc(end,1)<V_t
    disp('No spike has occurred! Try larger input current I_inj.');
    return
end
[t, q, ncond] = calc_iPRC(t,lc,param,I_inj);  % infinitesimal PRC
toc

%%% plot
figure;
subplot(3,2,1);
plot(t,lc(:,1),'linewidth',1.5);
title('limit cycle trajectory');
ylabel('V (mV)');
subplot(3,2,3);
plot(t,0.001*lc(:,2),'linewidth',1.5);
ylabel('w (nA)');
subplot(3,2,5);
plot([t(1) t(end)], [0 0], 'k--');
plot(t,pert_V*q(:,1),'linewidth',1.5) %iPRC multiplied by perturb. strength
title(['PRC for perturbations of strength ' num2str(pert_V) 'mV']);
ylabel('phase shift (ms)');
xlabel('t (ms)');

subplot(3,2,[2 6]); hold on  % phase plane plot
V_range = -70:0.1:-35;
V_nullcline = -g_l*(V_range-E_l)+g_l*delta_T*exp((V_range-V_t)/delta_T);
w_nullcline = a*(V_range-E_l);

% vector field
Vspan = -65:4:-40;
wspan = 300:10:400;
[Vf, wf] = meshgrid(Vspan,wspan);
DV = (-g_l*(Vf-E_l)+g_l*delta_T*exp((Vf-V_t)/delta_T)-wf+I_inj)/C;
Dw = (a*(Vf-E_l)-wf)/tau_w;

zidx = [];
for j=1:length(q(:,1))-2
    if q(j,1)<0 && q(j+1,1)>0
      zidx = [zidx j];
    elseif q(j,1)<0 && q(j+1,1)==0 && q(j+2,1)>0
      zidx = [zidx j+1];
    end
end

% nullclines and limit cycle
plot(V_range,V_nullcline+I_inj,'k','linewidth',1.5);
plot(V_range,w_nullcline,'k--','linewidth',1.5);
[~,midx] = max(q(:,1));
plot(lc(:,1),lc(:,2),'g','linewidth',2);
plot(lc(1,1),lc(1,2),'go','markerfacecolor','g'); 
plot(lc(midx,1),lc(midx,2),'gd','markerfacecolor','g'); 
plot(lc(zidx,1),lc(zidx,2),'gs','markerfacecolor','g');
legend('V-nullcline','w-nullcline','aEIF trajectory','reset',...
       'PRC maximum','PRC 0-crossing');

% plot vector field 
quiver(Vf,wf,DV,Dw,3,'color',[0.5 0.5 0.5]);
axis([Vspan(1) Vspan(end) wspan(1) wspan(end)]);
set(gca,'xtick',[Vspan(1) Vspan(end)],'ytick',[wspan(1) wspan(end)],...
    'fontsize',12,'box','on');
xlabel('V'); ylabel('w');
title('state space');
axis square