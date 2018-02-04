function [t, q, ncond] = calc_iPRC(tp,lc,param,I_inj)
% solves the adjoint system dq/dt = -Df(lc)' * q, with q := [q_V; q_w]  
% subject to the normalization condition q(t)'*f(lc(t)) = 1 and ensuring 
% q_w(0) = q_w(T), as described in Ladenbauer et al., PLOS Comput. Biol. 2012
% written by Josef Ladenbauer

C = param(1);   % membrane capacitance [pF]
g_l = param(2);   % leak conductance [nS]
E_l = param(3);   % leak reversal potential [mV]
V_t = param(4);   % spike threshold [mV]
delta_T = param(5);   % slope factor [mV]
tau_w = param(6);    % adaptation time constant [ms]
a = param(7);   % subthreshold adaptation [nS]
b = param(8);   % spike-triggered adaptation [pA]

ncond = zeros(length(tp),1);

if a==0 && b==0  
    
    z = f_eval(lc(end,:));
    q(length(tp)) = 1/z(1);   % initial point, normalized    
    for i=length(tp):-1:2   % backward integration
        J = (-g_l+g_l*exp((lc(i,1)-V_t)/delta_T))/C;
        q(i-1) = q(i) + (tp(i)-tp(i-1)) * J * q(i);   % Euler method
    end
    z = f_eval(lc(1,:));
    q = q / (q(1)*z(1));
    q = [q' zeros(length(q),1)];
    for i=1:length(tp) 
        ncond(i) = q(i,:)*f_eval(lc(i,:));
    end
    t = tp;

else
    
    z = f_eval(lc(end,:));

    q(length(tp),1) = z(1);  % initial point
    q(length(tp),2) = z(2);

    opts = optimset('TolFun',1e-7);
    q2_T = fzero(@ivp,-0.5,opts);
    % normalized initial point
    q(end,1) = (1 - q2_T*z(2)) / z(1);   
    q(end,2) = q2_T;
    
    for i=1:length(tp)
        ncond(i) = q(i,:)*f_eval(lc(i,:));
    end
    t = tp;

end


function q2_diff = ivp(q2_T)
  q(length(tp),1) = (1 - q2_T*z(2)) / z(1);   % normalized initial point
  q(length(tp),2) = q2_T;
  for j=length(tp):-1:2   % backward integration       
    J = [(-g_l+g_l*exp((lc(j,1)-V_t)/delta_T))/C -1/C; ...
         a/tau_w -1/tau_w];   % Jacobian at lc point: Df(lc(tp(i)))
    q(j-1,:) = q(j,:) + (tp(j)-tp(j-1)) * q(j,:) * J;   % Euler method      
  end
  q2_diff = q(end,2) - q(1,2);
end

function z = f_eval(x)
  z = [(-g_l*(x(1)-E_l)+g_l*delta_T*exp((x(1)-V_t)/delta_T)-x(2)+I_inj)/C; 
      (a*(x(1)-E_l)-x(2))/tau_w];
end
  
end