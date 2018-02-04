% script for computing the limit cycle trajectory and infinitesimal phase 
% response curve (iPRC) for an adaptive exponential integrate-and-fire 
% neuron, as well as the interaction function H between two such neurons 
% coupled through a conductance-based synapse model;
% written by Josef Ladenbauer

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
params = [C g_l E_l V_t delta_T tau_w a b V_r];
[t, lc] = aEIF_lc(params,I_inj,V_r,w_init,T_max);  % limit cycle
if lc(end,1)<V_t
    disp('No spike has occurred! Try larger input current I_inj.');
    return
end
[t, q, ncond] = calc_iPRC(t,lc,params,I_inj);  % infinitesimal PRC
toc

%%% plot
figure;
subplot(4,2,1);
plot(t,lc(:,1),'linewidth',1.5);
title('limit cycle trajectory: V(t)');
ylabel('V (mV)');
subplot(4,2,3);
plot(t,0.001*lc(:,2),'linewidth',1.5);
title('limit cycle trajectory: w(t)');
ylabel('w (nA)');
subplot(4,2,5);
plot([t(1) t(end)], [0 0], 'k--');
plot(t,pert_V*q(:,1),'linewidth',1.5) %iPRC multiplied by perturb. strength
title(['PRC for perturbations of strength ' num2str(pert_V) 'mV']);
ylabel('phase shift (ms)');
% xlabel('t (ms)');

subplot(4,2,[2 6]); hold on  % phase plane plot
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

% compute interaction function H for two identical aEIF neurons connected
% by a conductance-based bi-exponential synapse model, using the iPRC q

lc_V = lc(:,1);
q_V = q(:,1);
tau_r = 0.1;  % rise time constant (ms)
tau_d = 1;  % decay time constant (ms)
delay = 10;  % propagation delay (ms)
E_syn = 0;  % synaptic reversal potential (mV)
g_syn = 1;  % coupling strength (nS)
phi = 0;  % phase shift between the neurons (ms)
syn_params = [tau_r tau_d delay E_syn g_syn/C];
[tp,H,pert,H_pt] = H_function(t,lc_V,q_V,t(end),syn_params,phi); 

subplot(4,2,7);
plot(tp,H,'k','linewidth',1.5); hold on
plot(phi,H_pt,'ko');
title('interaction function H for 2 coupled neurons');
xlabel('t (ms)');