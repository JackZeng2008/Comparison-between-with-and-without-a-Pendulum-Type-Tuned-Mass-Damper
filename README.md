# Comparison-between-with-and-without-a-Pendulum-Type-Tuned-Mass-Damper
It is Matlab code that examines a pendulum-type mass damper's effectiveness of mitigating displacement and acceleration of a building.
clear; clc;

% Load earthquake data
fname = 'CO2_RSN6_IMPVALL.I_I-ELC180.txt';
data = load(fname);
tg = data(:,1); % Time vector
ag = data(:,2); % Acceleration data (in g)
scaleFactor = 9.81; % Conversion factor from g to m/s²
ag_ms2 = ag .* scaleFactor; % Convert to m/s²
agFunc = @(t) interp1(tg, ag_ms2, t, 'linear', 'extrap'); % Interpolation function

% System parameters
M = 2.0e7;      % Main structure mass (kg)
K = 3.16e9;     % Main structure stiffness (N/m)
m_tmd = 8.0e5;  % TMD mass (kg)
g = 9.81;       % Gravitational acceleration (m/s²)
zeta = 0.02;    % Fixed damping ratio
C = 2 * zeta * sqrt(K * M); % Structural damping coefficient

% Optimize TMD parameters ================================================
% Calculate natural properties of main structure
omega_n = sqrt(K/M);     % Natural circular frequency (rad/s)
f_n = omega_n/(2*pi);    % Natural frequency (Hz)
T_n = 1/f_n;             % Natural period (s)

% Optimize TMD frequency (slightly below structure frequency)
tuning_ratio = 0.98;     % Tuning ratio (0.95-1.0)
omega_tmd = omega_n * tuning_ratio;

% Optimal pendulum length
l = g / omega_tmd^2;     % Original 5m, optimized to ~3.97m

% Optimal damping ratio (Den Hartog formula)
mu = m_tmd / M;          % Mass ratio (10%)
zeta_tmd_opt = sqrt(3*mu/(8*(1+mu))); % ~7.5%

% Optimal damping coefficient (convert to torque coefficient)
cd_linear = 2 * zeta_tmd_opt * (m_tmd * omega_tmd); % Linear damping
cd = cd_linear * l^2;    % Original 1000, optimized to ~1260
% ==========================================================

% Time settings
tspan = [tg(1) tg(end)]; % Simulation time span
x0 = [0; 0; 0; 0]; % Initial conditions [displacement; velocity; angle; angular velocity]

% Solve without TMD
odefun_noTMD = @(t,x) noTMD_ode(t, x, M, C, K, agFunc);
opts = odeset('RelTol',1e-6,'AbsTol',1e-9); % ODE solver options
[t_noTMD, x_noTMD] = ode45(odefun_noTMD, tspan, [0;0], opts);

% Calculate acceleration without TMD
numSteps = numel(t_noTMD);
acc_rel_noTMD = zeros(numSteps,1); % Relative acceleration
for ii = 1:numSteps
    dxdt = odefun_noTMD(t_noTMD(ii), x_noTMD(ii,:)');
    acc_rel_noTMD(ii) = dxdt(2); % Second state is velocity
end
ag_vect = agFunc(t_noTMD); % Ground acceleration
acc_abs_noTMD = acc_rel_noTMD + ag_vect; % Absolute acceleration

% Solve with TMD
odefun_TMD = @(t,x) TMD_ode(t, x, M, C, K, m_tmd, l, cd, g, agFunc);
[t_TMD, x_TMD] = ode45(odefun_TMD, tspan, x0, opts);

% Calculate acceleration with TMD
numSteps = numel(t_TMD);
acc_rel_TMD = zeros(numSteps,1);
for ii = 1:numSteps
    dxdt = odefun_TMD(t_TMD(ii), x_TMD(ii,:)');
    acc_rel_TMD(ii) = dxdt(2);
end
ag_vect = agFunc(t_TMD);
acc_abs_TMD = acc_rel_TMD + ag_vect;

% Extract structural displacements
disp_noTMD = x_noTMD(:,1); % Displacement without TMD
disp_TMD = x_TMD(:,1);     % Displacement with TMD

% Plot: Displacement comparison
figure('Name','Structural Displacement: With vs Without TMD');
hold on;
plot(t_noTMD, disp_noTMD, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Without TMD');
plot(t_TMD, disp_TMD, 'r-', 'LineWidth', 1.5, 'DisplayName', 'With TMD');
xlabel('Time (s)');
ylabel('Displacement (m)');
title(['Displacement Comparison (\zeta_{str}=0.02, \zeta_{TMD}=', num2str(zeta_tmd_opt*100,2),'%)']);
legend('show', 'Location', 'best');
grid on;
hold off;

% Plot: Acceleration comparison (new section)
figure('Name','Structural Acceleration: With vs Without TMD');
plot(t_noTMD, acc_abs_noTMD, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Without TMD');
hold on;
plot(t_TMD, acc_abs_TMD, 'r-', 'LineWidth', 1.5, 'DisplayName', 'With TMD');
xlabel('Time (s)');
ylabel('Absolute Acceleration (m/s^2)');
title(['Acceleration Comparison (\zeta_{str}=0.02, \zeta_{TMD}=', num2str(zeta_tmd_opt*100,2),'%)']);
legend('Location', 'best');
grid on;
hold off;

% Calculate peak responses
peakDisp_noTMD = max(abs(disp_noTMD));
peakAccel_noTMD = max(abs(acc_abs_noTMD));
peakDisp_TMD = max(abs(disp_TMD));
peakAccel_TMD = max(abs(acc_abs_TMD));

% Print optimized parameters
fprintf('===== TMD Optimization Parameters =====\n');
fprintf('Structure natural period: %.3f s\n', T_n);
fprintf('TMD tuning ratio: %.3f\n', tuning_ratio);
fprintf('Optimal pendulum length: %.3f m\n', l);
fprintf('Optimal TMD damping ratio: %.1f%%\n', zeta_tmd_opt*100);
fprintf('Optimal damping coefficient: %.1f N·m·s/rad\n', cd);
fprintf('======================================\n');

% Print results
fprintf('\n===== With vs Without TMD Comparison =====\n');
fprintf('                     Without TMD        With TMD       Reduction\n');
fprintf('Peak Disp (m)     %14.6e    %14.6e    %6.2f%%\n', ...
        peakDisp_noTMD, peakDisp_TMD, 100*(peakDisp_noTMD - peakDisp_TMD)/peakDisp_noTMD);
fprintf('Peak Accel (m/s²) %14.6e    %14.6e    %6.2f%%\n', ...
        peakAccel_noTMD, peakAccel_TMD, 100*(peakAccel_noTMD - peakAccel_TMD)/peakAccel_noTMD);
fprintf('===================================================\n');

% ODE function for system without TMD
function dxdt = noTMD_ode(t, x, M, C, K, agFunc)
    dxdt = zeros(2,1);
    ag = agFunc(t); % Get ground acceleration at time t
    dxdt(1) = x(2); % Velocity
    dxdt(2) = (-C*x(2) - K*x(1) - M*ag)/M; % Acceleration
end

% ODE function for system with TMD
function dxdt = TMD_ode(t, x, M, C, K, m, l, cd, g, agFunc)
    dxdt = zeros(4,1);
    ag = agFunc(t); % Ground acceleration
    
  % Mass matrix
    massMatrix = [M,      m*l;  % Structure mass and coupling term
                 m*l, m*l^2];   % TMD inertia terms
    
  % Force vector
    forceVector = [-C*x(2) - K*x(1) - M*ag;  % Structure forces
                   -cd*x(4) - m*g*l*x(3) - m*l*ag]; % TMD forces
    
  % Solve for accelerations
    accelerations = massMatrix \ forceVector;
    
  % State derivatives
    dxdt(1) = x(2);             % Structure velocity
    dxdt(2) = accelerations(1); % Structure acceleration
    dxdt(3) = x(4);             % Pendulum angular velocity
    dxdt(4) = accelerations(2); % Pendulum angular acceleration
end
