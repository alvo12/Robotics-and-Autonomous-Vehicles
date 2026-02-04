clc; clear; close all;

% Vehicle Parameters
m  = 1.0;     
lf = 0.075;
lr = 0.075;
Iz = 0.004;
Cf = 2000;
Cr = 2000;

% Speeds
V = [4 8 12];  % m/s

% LQR weights
Q = diag([10 10 1 1]);
R = 1;

% Time settings
dt = 0.004;
t  = 0:dt:5;

% Initial condition
x0 = [0.05; 5*pi/180; 0; 0];

% Prepare figures
figure(1); hold on; grid on; title('Lateral Position y(t)');
xlabel('Time (s)'); ylabel('y (m)');

figure(2); hold on; grid on; title('Heading Error \psi(t)');
xlabel('Time (s)'); ylabel('\psi (rad)');

figure(3); hold on; grid on; title('Lateral Velocity v_y(t)');
xlabel('Time (s)'); ylabel('v_y (m/s)');

figure(4); hold on; grid on; title('Yaw Rate r(t)');
xlabel('Time (s)'); ylabel('r (rad/s)');

% Loop over speeds
for k = 1:length(V)
    
    v0 = V(k);
    
    % Lateral dynamics
    A_lat = [-(Cf+Cr)/(m*v0), (lr*Cr-lf*Cf)/(m*v0)-v0;
              (lr*Cr-lf*Cf)/(Iz*v0), -(lf^2*Cf+lr^2*Cr)/(Iz*v0)];
    
    B_lat = [Cf/m; lf*Cf/Iz];
    
    % Augmented model
    A_e = [0 v0 1 0;
           0 0  0 1;
           0 0 A_lat(1,1) A_lat(1,2);
           0 0 A_lat(2,1) A_lat(2,2)];
    
    B_e = [0;0;B_lat];
    
    % LQR gain
    K = lqr(A_e, B_e, Q, R);
    
    % Closed-loop system
    A_cl = A_e - B_e*K;
    sys = ss(A_cl, [], eye(4), []);
    
    % Simulate
    [~, ~, x] = initial(sys, x0, t);
    
    % Plot states
    figure(1); plot(t, x(:,1), 'LineWidth', 2);
    figure(2); plot(t, x(:,2), 'LineWidth', 2);
    figure(3); plot(t, x(:,3), 'LineWidth', 2);
    figure(4); plot(t, x(:,4), 'LineWidth', 2);
    
end

% Legends
figure(1); legend('v = 4 m/s','v = 8 m/s','v = 12 m/s');
figure(2); legend('v = 4 m/s','v = 8 m/s','v = 12 m/s');
figure(3); legend('v = 4 m/s','v = 8 m/s','v = 12 m/s');
figure(4); legend('v = 4 m/s','v = 8 m/s','v = 12 m/s');
