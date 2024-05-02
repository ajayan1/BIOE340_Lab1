
% Parameters
k_plus = 5e-4; % rate that O2 binds to Hb
k_minus = 1.8e5*k_plus; % rate that O2 leaves Hb
sigma_O2_B = 1.2e-6; %PO2 in the blood


% Initial conditions
Thb = 2.2e-3; %total concentration of Hb (M)
HB4 = 1.7e-3;
HBm=Thb-HB4;
PO2_B = 40;

% Time span
tspan = 0:0.1:2;

% Solve differential equations
[t, y] = ode45(@(t, y) gas_exchange_eqns(t, y, k_plus, k_minus, sigma_O2_B,Thb,HB4,PO2_B,HBm), tspan, [HBm; HB4; PO2_B]);

% Plot results
figure;
subplot(3, 1, 1);
plot(t, y(:, 1), 'b', 'LineWidth', 2);
xlabel('Time');
ylabel('HBm');
title('HBm vs. Time');

subplot(3, 1, 2);
plot(t, y(:, 2), 'r', 'LineWidth', 2);
xlabel('Time');
ylabel('HB4');
title('HB4 vs. Time');

subplot(3, 1, 3);
plot(t, y(:, 3), 'g', 'LineWidth', 2);
xlabel('Time');
ylabel('PO2_B');
title('PO2_B vs. Time');



function dydt = gas_exchange_eqns(t, y, k_plus, k_minus, sigma_O2_B, Thb, HB4,PO2_B,HBm)
    
    m=3.6; % number of O2 molecules needed to bidn to Hb to form Hb4
    % Extract variables
    HBm = y(1);
    HB4 = y(2);
    PO2_B = y(3);

    % Differential equations
    dHB4dt = k_plus*(Thb-HB4)*PO2_B - k_minus*HB4;
    dHBmdt = -dHB4dt;
    dPO2_Bdt = m* (-k_plus*(Thb-HB4)*PO2_B + k_minus*HB4)/sigma_O2_B;

    % Return derivatives
    dydt = [dHBmdt; dHB4dt; dPO2_Bdt];

end
