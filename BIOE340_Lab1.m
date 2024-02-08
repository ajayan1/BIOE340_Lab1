clc;close all;clear;
% variables 
C = 1;         % Membrane capacitance
Ena = 50;      % Sodium reversal potential
Ek = -77;      % Potassium reversal potential
El = -49;      % Leak reversal potential
Gna = 120;     % Sodium conductance 
Gk = 36;       % Potassium conductance 
Gl = 0.3;      % Leak conductance 
    
%initial n, m, h, V
V0 = -65;      % Resting membrane potential
n0 = 0.337;    % Initial value for n
m0 = 0.061;    % Initial value for m
h0 = 0.552;    % Initial value for h
y0 = [V0, n0, h0, h0];
tspan = [0 6];
%ODE45 
f = @(t,y) HodgHux(t,y);
options = odeset('MaxStep',0.01);
[t,y] = ode45(f,tspan,y0,options);

% plots 
 figure;
    
% Plot membrane potential over time
subplot(2, 1, 1);
plot(t, y(:, 1), 'b-');
title('Membrane Potential over Time');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
grid on;
    
% Plot gating parameters (n, m, h) over time
subplot(2, 1, 2);
plot(t, y(:, 2), 'r-','DisplayName', 'n');
hold on;
plot(t, y(:, 3), 'k-','DisplayName', 'm');
plot(t, y(:, 4), 'm-', 'DisplayName', 'h');
title('Gating Parameters over Time');
xlabel('Time (ms)');
ylabel('Gating Parameter');
legend('Location', 'Best');
grid on;
%%
function lab2 = HodgHux(t,y)
% State the variables you will be using 
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);
    C = 1;         % Membrane capacitance
    Ena = 50;      % Sodium reversal potential
    Ek = -77;      % Potassium reversal potential
    El = -49;      % Leak reversal potential
    Gna = 120;     % Sodium conductance 
    Gk = 36;       % Potassium conductance 
    Gl = 0.3;      % Leak conductance 
    vm = V+62;
%write out the functions for: alpha h, beta h, beta n, and beta m
    alpha_h = 0.07 .* exp(-vm ./ 20);
    beta_h = (exp((30-vm)/10) + 1).^-1;
    alpha_n = 0.01 .* (10-vm) ./ (exp((10-vm)./10)-1);
    beta_n = 0.125 .* exp(-vm ./ 80);
    alpha_m = 0.1.*(25-vm) ./ (exp((25-vm)./10)-1);
    beta_m = 4 .* exp(-vm./18);
%loop for alpha n if Vd == 10 and alpha m if Vd == 25 
    if V == 25
        alpha_m = 1;
    end
    if V ==10
        alpha_n = 0.1;
    end
%write out the ODE for n, m, and h
    dn = alpha_n.*vm .* (1 - n) - beta_n .* n;
    dm = alpha_m .* (1 - m) - beta_m .* m;
    dh = alpha_h .* (1 - h) - beta_h .* h;
%write out the conductance equations 
    Gk = vm.*n^4;
    Gna = vm .* m.^3.*h;
%write out the ODE for V
    dV = (-Gna.*(V-Ena)-Gk.*(V-Ek)-Gl.*(V-El))./C;
% loop for dV @ t>1 && t<1.1
    if t>1 && t<1.1
        dV = 430;
    end
lab2 = [dV;dn;dm;dh];
end
