function [tp, H, pert_save, H_pt] = H_function(t,lc,q,T,syn_params,phi)
% calculates the interaction function H as a function of phase difference
% for a given limit cycle lc (V component) and iPRC q (V component) of two
% identical neurons, using a conductance-based bi-exponential synapse model 

tau_r = syn_params(1);
tau_d = syn_params(2);
delay = syn_params(3);
E_syn = syn_params(4);
cs = syn_params(5);  % coupling strength (nS/pF)

%%% calculate interaction function by applying phase reduction theory

% interpolate lc and q on uniform time grid
% (need uniform time grid to obtain odd part of H below)
res = 500;
dt = T/res;
tp = 0:dt:T; 
lc_int = interp1(t,lc,tp,'linear');
q_int = interp1(t,q,tp,'linear');
[~,idx] = min(abs(tp-phi));
H = zeros(1,length(tp));

for i=1:length(tp)
    % vary the phase of neuron 2: theta_1 (=tp) + phi - delay
    % where phi = theta_2 - theta_1 (phase difference) given by tp(i)
    th_j = mod(tp+tp(i)-delay,T);
    
    % conductance based synapse I_syn(t) = g_syn *s(t) * (E_syn-V(t)) with
    % s(t) given by bi-exponential function, peak normalized to 1:
    s = exp(-th_j/tau_d) - exp(-th_j/tau_r);
    s = s/max(s);
    pert = cs * s .* (E_syn - lc_int);
    
    if i==idx
      pert_save = pert;
    end
    % average over the oscillation cycle to obtain interaction functions H
    H(i) = mean(q_int.*pert);  % H of neuron receiving syn. input
end

H_pt = H(idx);

end