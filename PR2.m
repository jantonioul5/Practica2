clc; clear; close all;

% Parámetros del sistema
sys.Ip = 0.0079;    % kg*m^2  Momento de inercia del péndulo
sys.Mc = 0.7031;    % kg      Masa del carrito
sys.lp = 0.3302;    % m       Distancia pivote-CG
sys.Mp = 0.23;      % kg      Masa del péndulo
sys.Fc = 0;       % N       Fuerza aplicada
sys.Bc = 4.3;    % Ns/m    Amortiguamiento carrito
sys.Bp = 0.0024; % Nms/rad Amortiguamiento péndulo
sys.g = 9.81;        % m/s^2   Gravedad

% Condiciones iniciales
t = deg2rad(1);    % Ángulo inicial rad
xc = 0;            % Posición inicial carro
t0 = 0;         % Velocidad angular inicial
xd = 0;             % Velocidad inicial carro

X0 = [xc; xd; t; t0];

% Tiempo de simulación
Tspan = [0 20];

% Resolver sistema no lineal
[T, Y] = ode45(@(t,X) pendulumDynamics(t,X,sys), Tspan, X0);

function dXdt = pendulumDynamics(~, X, sys)
    % Variables de estado
    x = X(1);          
    x_dot = X(2);
    theta = X(3);     
    theta_dot = X(4);

    % Parámetros
    Ip = sys.Ip; Mc = sys.Mc; Mp = sys.Mp; lp = sys.lp;
    Fc = sys.Fc; Bc = sys.Bc; Bp = sys.Bp; g = sys.g;

    % Denominador común
    D = (Mc+Mp)*Ip + Mc*Mp*lp^2 + Mp^2*lp^2*sin(theta)^2;

    % Aceleración del carro
    xdd = ((Ip+Mp*lp^2)*Fc + Mp^2*lp^2*g*sin(theta)*cos(theta) ...
             - (Ip+Mp*lp^2)*Bc*x_dot - (Ip*Mp*lp - Mp^2*lp^3)*theta_dot^2*sin(theta) ...
             - Mp*lp*Bp*theta_dot*cos(theta)) / D;

    % Aceleración angular del péndulo
    tdd = ((Mc+Mp)*Mp*g*lp*sin(theta) - (Mc+Mp)*Bp*theta_dot ...
                 + Fc*Mp*lp*cos(theta) - Mp^2*lp^2*theta_dot^2*sin(theta)*cos(theta) ...
                 - Bc*Mp*lp*x_dot*cos(theta)) / D;

    dXdt = [x_dot; xdd; theta_dot; tdd];
end

% Graficación

figure('Color',[0.95 0.95 0.95],'Position',[100 100 900 600]);

subplot(2,2,1);
plot(T,Y(:,1),'Color',[0 0.45 0.74],'LineWidth',1.8);
xlabel('Tiempo [s]'); ylabel('x [m]');
title('Posición del carro'); grid on;

subplot(2,2,3);
plot(T,Y(:,2),'Color',[0.47 0.67 0.19],'LineWidth',1.8);
xlabel('Tiempo [s]'); ylabel('v [m/s]');
title('Velocidad del carrito'); grid on;

subplot(2,2,2);
plot(T,Y(:,3),'Color',[0.85 0.33 0.10],'LineWidth',1.8);
xlabel('Tiempo [s]'); ylabel('\theta [rad]');
title('Ángulo del péndulo'); grid on;

subplot(2,2,4);
plot(T,Y(:,4),'Color',[0.49 0.18 0.56],'LineWidth',1.8);
xlabel('Tiempo [s]'); ylabel('\omega [rad/s]');
title('Velocidad angular del péndulo'); grid on;

% Linealización para espacio de estados

Mc = sys.Mc; Mp = sys.Mp; Ip = sys.Ip; lp = sys.lp;
g = sys.g; Bc = sys.Bc; Bp = sys.Bp;

A = [0 1 0 0;
     0 -Bc/Mc (Mp*g)/Mc 0;
     0 0 0 1;
     0 -Bc*Mp*lp/(Mc*Ip) g*(Mc+Mp)*lp/Ip -Bp/Ip];

B = [0; 1/Mc; 0; lp/Ip];
C = eye(4);
D = zeros(4,1);

disp('Matriz A del sistema linealizado:'); disp(A)
disp('Matriz B:'); disp(B)