// 04/16/2014 - ENGR 2450 - Meine, Joel
// Lorenz Equations (Pg. 816 | Chapra & Canale)

function [Z] = F(t,Y)
    Z(1) = -o*Y(1) + o*Y(2);
    Z(2) = r*Y(1) - Y(2) - Y(1)*Y(3);
    Z(3) = -b*Y(3) + Y(1)*Y(2);
endfunction
o = 10; b = 2.666667; r = 28;
dt = 0.01; // 0.5 for table values; 0.01 for plots
t = [0.0:dt:20.0];
t0 = 0;
Y0 = [5;5;5];
Y = ode("rk",Y0,t0,t,F);
disp("t   x   y   z");
disp([t' Y']);
y1 = Y(1,:); y2 = Y(2,:); y3 = Y(3,:);
scf(); plot(t,y1,'b-',t,y2,'r-',t,y3,'g-');
legend('x','y','z',1);
xtitle('Signal Plots','x','y','z');
scf(); plot(y1,y2);
xtitle('Phase Portrait','x','y');
scf(); plot(y2,y3);
xtitle('Phase Portrait','y','z');
scf(); plot(y3,y1);
xtitle('Phase Portrait','z','x');
