// 04/16/2014 - ENGR 2450 - Meine, Joel
// Problem 25.27

function [dydt] = f(t,y)
    dydt = -k*sqrt(y);
endfunction;
t0 = 0; y0 = 3; tmax = 57.5; h = 0.5;
k = 0.06;
t = [t0:h:tmax];
y = ode("rk",y0,t0,t,f);
disp("t   y");
disp("--------");
disp([t' y']);
function dydt = fe(t)
    dydt = ((k^2)/4)*((2/k)*sqrt(y0)-t)^2;
endfunction;
te = [t0:h/10:tmax];
n = length(te);
for i = 1:n
    ye(i)=fe(te(i));
end;
plot(t,y,'ro',te,ye,'-b'); legend('RK4','Exact',3);
xtitle('Solving dy/dt = -k*sqrt(y)','t','y');
