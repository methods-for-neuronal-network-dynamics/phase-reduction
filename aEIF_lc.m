function [t, lc] = aEIF_lcf(param,I_inj,V_init,w_init,Tmax)
%%% Computes the stable limit cycle of an oscillating aEIF model neuron
%%% using the parameters from the vector 
%%% param = [C  g_l  E_l  V_t  delta_T  tau_w  a  b  V_r  I_inj].
%%% Returns limit cycle (t, lc) where lc = [V;w] at time points tp in [0 T]

%% parameters
C = param(1);   % membrane capacitance [pF]
g_l = param(2);   % leak conductance [nS]
E_l = param(3);   % leak reversal potential [mV]
V_t = param(4);   % spike threshold [mV]
delta_T = param(5);   % slope factor [mV]
tau_w = param(6);    % adaptation time constant [ms]
a = param(7);   % subthreshold adaptation [nS]
b = param(8);   % spike-triggered adaptation [pA]
V_r = param(9);   % reset membrane voltage [mV] 
V_cut = -30;   % cutoff voltage [mV], somewhere between V_t and 0

rel_w_dev_max = 1e-4;
max_iter = 50;

lc0 = [V_init; w_init];

options = odeset('Events',@spike,'InitialStep',0.05,'MaxStep',0.05);
count = 1;  rel_w_dev = 1;
while count<max_iter && rel_w_dev>rel_w_dev_max
    [t,lc,~,~,~] = ode15s(@odesys,[0 Tmax],lc0,options);
    if a>0 || b>0
        rel_w_dev = abs((lc(1,2)-lc(end,2)-b)/lc(1,2));
    else
        rel_w_dev = 0;
    end
    lc0 = [V_r; lc(end,2)+b];
    count = count+1;
end

if rel_w_dev>rel_w_dev_max
    disp(['Inaccurate limit cycle ' ...
          '==> try again with increased max_iter in aEIF_lcf.m'])
    return
end


%% ODE system
function dlc = odesys(t,lc)
  dlc = [0;0];
  dlc(1) = (-g_l*(lc(1)-E_l) + g_l*delta_T*exp((lc(1)-V_t)/delta_T)- ...
           lc(2) + I_inj)/C;
  dlc(2) = (a*(lc(1)-E_l)-lc(2))/tau_w;
end

%% spike events
function [lookfor, stop, direction] = spike(t,lc)
  lookfor = lc(1)-V_cut;   % stop when V = cutoff voltage
  stop = 1; 
  direction = 1;   % cutoff voltage reached from below
end
end