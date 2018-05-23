// 04/16/2014 - ENGR 2450 - Meine, Joel
// Problem 28.47

function [Z] = F(t,Y)
    Z(1) = Y(2);
    Z(2) = (Fo*sin(w*t)-a*abs(Y(2))*Y(2)-k*Y(1))/m;
endfunction
m = 2; a = 5; k = 6; Fo = 2.5; w = 0.5;
dt = 0.25;
t = [0.0:dt:15.0];
t0 = 0;
Y0 = [1;0];
Y = ode("rk",Y0,t0,t,F);
disp("t   x   dx/dt");
disp([t' Y']);
y1 = Y(1,:); y2 = Y(2,:);
scf(); plot(t,y1,'b-',t,y2,'r-');
legend('x','dx/dt',1);
xtitle('Signal Plots','x','dx/dt');
scf(); plot(y1,y2);
xtitle('Phase Portrait','x','dx/dt');
