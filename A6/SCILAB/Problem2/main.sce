// 04/16/2014 - ENGR 2450 - Meine, Joel
// Problem 25.21

function [dpdt] = f(t,p)
    dpdt = kgm*(1-(p/pmax))*p;
endfunction;
t0 = 1950; p0 = 2555; tn = 2000; h = 10;
kgm = 0.026; pmax = 12000;
t = [t0:h:tn];
p = ode("rk",p0,t0,t,f);
disp("t   p");
disp("--------");
disp([t' p']);
function dpdt = fe(t)
    dpdt = pmax / (1-(1-(pmax/p0))*exp(-kgm*(t-t0)));
endfunction;
te = [t0:h/10:tn];
n = length(te);
for i = 1:n
    pe(i)=fe(te(i));
end;
plot(t,p,'ro',te,pe,'-b'); legend('RK4','Exact',3);
xtitle('Solving dp/dt = kgm*(1-(p/pmax))*p','t','p');
