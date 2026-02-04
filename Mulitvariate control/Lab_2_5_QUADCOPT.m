clc; clear; close all;

%% ===================== QUADCOPTER PARAMETERS ============================

% Moments of inertia (kg·m^2)
Ix = 0.005;
Iy = 0.005;
Iz = 0.009;

%% ===================== LINEARIZED MODEL ================================

A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0];

B = [0     0     0;
     0     0     0;
     0     0     0;
     1/Ix  0     0;
     0     1/Iy  0;
     0     0   1/Iz];

%% ===================== LQR DESIGN =====================================

Q = diag([100 100 100 10 10 10]);   % State penalties
R = diag([1 1 1]);                  % Control penalties

K = lqr(A, B, Q, R);                % LQR gain

disp('LQR Gain Matrix:');
disp(K)

%% ===================== SIMULATION SETTINGS =============================

dt = 0.002;          % Time step (s)
T  = 5;              % Total simulation time (s)
t  = 0:dt:T;
N  = length(t);

%% ===================== INITIAL CONDITIONS ==============================

% Initial angles: 10 deg roll, pitch, yaw
x = zeros(6, N);
x(:,1) = [10*pi/180;
          10*pi/180;
          10*pi/180;
          0;
          0;
          0];

%% ===================== RK4 SIMULATION LOOP =============================

for k = 1:N-1
    
    % LQR control law
    u = -K * x(:,k);
    
    % RK4 integration
    k1 = quad_dynamics(x(:,k), u, Ix, Iy, Iz);
    k2 = quad_dynamics(x(:,k) + dt/2 * k1, u, Ix, Iy, Iz);
    k3 = quad_dynamics(x(:,k) + dt/2 * k2, u, Ix, Iy, Iz);
    k4 = quad_dynamics(x(:,k) + dt   * k3, u, Ix, Iy, Iz);
    
    x(:,k+1) = x(:,k) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

%% ===================== PLOTTING ========================================

figure('Color','w')

subplot(3,1,1)
plot(t, x(1,:)*180/pi, 'LineWidth', 2)
grid on
ylabel('\phi (deg)')
title('LQR Stabilized Quadcopter Attitude')

subplot(3,1,2)
plot(t, x(2,:)*180/pi, 'LineWidth', 2)
grid on
ylabel('\theta (deg)')

subplot(3,1,3)
plot(t, x(3,:)*180/pi, 'LineWidth', 2)
grid on
ylabel('\psi (deg)')
xlabel('Time (s)')

%% ===================== ANGULAR RATES ===================================

figure('Color','w')

subplot(3,1,1)
plot(t, x(4,:), 'LineWidth', 2)
grid on
ylabel('p (rad/s)')
title('Angular Rates')

subplot(3,1,2)
plot(t, x(5,:), 'LineWidth', 2)
grid on
ylabel('q (rad/s)')

subplot(3,1,3)
plot(t, x(6,:), 'LineWidth', 2)
grid on
ylabel('r (rad/s)')
xlabel('Time (s)')

%% ===================== NONLINEAR DYNAMICS FUNCTION ======================

function dx = quad_dynamics(x, u, Ix, Iy, Iz)

phi   = x(1);
theta = x(2);
psi   = x(3);
p     = x(4);
q     = x(5);
r     = x(6);

tau_phi   = u(1);
tau_theta = u(2);
tau_psi   = u(3);

dx = zeros(6,1);

% Kinematics
dx(1) = p;
dx(2) = q;
dx(3) = r;

% Dynamics (Euler rigid body equations)
dx(4) = ((Iy - Iz)/Ix) * q * r + tau_phi   / Ix;
dx(5) = ((Iz - Ix)/Iy) * p * r + tau_theta / Iy;
dx(6) = ((Ix - Iy)/Iz) * p * q + tau_psi   / Iz;

end
