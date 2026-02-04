X0 = [deg2rad(1); 0; 0; 0]; 
T = [0 10]; 
[t, x] = ode45(@Dinamica, T, X0);

figure(1) 
plot(t, x); 
grid on;
title('Espacio de estados');
legend('\alpha', '\alpha_p', 'x_c', 'x_p');

function dx = Dinamica(t, x)
 
    Ip = 0.0079; 
    Mc = 0.7031; 
    lp = 0.3302; 
    Mp = 0.23;
    Beq = 4.3; 
    Bp = 0.0024; 
    g = 9.81; 
    Fc = 0;

    dx = zeros(4,1); 
    
    den = (Mc + Mp)*Ip + Mc*Mp*lp^2 + Mp^2*lp^2*sin(x(1))^2;

    % x(1)=alpha, x(2)=d_alpha, x(3)=xc, x(4)=d_xc
    dx(1) = x(2);
    
    dx(2) = ((Mc + Mp)*Mp*g*lp*sin(x(1)) - (Mc + Mp)*Bp*x(2) + Fc*Mp*lp*cos(x(1)) ...
            - Mp^2*lp^2*x(2)^2*sin(x(1))*cos(x(1)) - Beq*Mp*lp*x(4)*cos(x(1))) / den;
    
    dx(3) = x(4); 
    
    dx(4) = ((Ip + Mp*lp^2)*Fc + Mp^2*lp^2*g*cos(x(1))*sin(x(1)) ...
            - (Ip + Mp*lp^2)*Beq*x(4) - (Ip*Mp*lp - Mp^2*lp^3)*x(2)^2*sin(x(1)) ...
            - Mp*lp*x(2)*cos(x(1))*Bp) / den;
end
