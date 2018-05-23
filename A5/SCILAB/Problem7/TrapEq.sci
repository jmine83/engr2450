function [T]=TrapEq(n, a, b, f)
    //This function implements the pseudocode for the
    //trapezoidal rule based on an equation y = f(x) as
    //described by Figure 22.1(a) in the Chapra-Canale
    //textbook, 6th edition
    h = (b-a)/n;
    x=a;
    suma = f(x);
    for i = 1:n-1
        x = x+h;
        suma = suma + 2*f(x);
    end;
    suma = suma + f(b);
    T = (b-a)*suma/(2*n);
endfunction;
