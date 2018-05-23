diary("Problems22_2And22_3.txt")
//Example of a script file used to run function Romberg
//Before running this program make sure that functions
//TrapEq and Romberg are loaded
//Problem22_2 --
a = 1; b = 2; maxiter = 50; ea = 0.5; //ea in percent
// define function to integrate as an inline function
function [y]=f222(x)
    y = (2*x+3/x)^2;
endfunction;
//Integrate using Romberg:
[I1,n1,iter1,ea1] = Romberg(a,b,maxiter,ea,f222)
//Problem22_3 --
a = 0; b = 2; maxiter = 50; ea = 0.5; //ea in percent
// define function to integrate as an inline function
function [z]=f223(t)
    z = exp(t)*sin(t)/(1+t^2);
endfunction;
//Integrate using Romberg:
[I2,n2,iter2,ea2] = Romberg(a,b,maxiter,ea,f223)
//End of script
diary(0)
