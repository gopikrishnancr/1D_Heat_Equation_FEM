function out = int_bool(f,a,b,n,t)
xval = linspace(a,b,5);    
bval = linspace(0,1,5);  % basis function values
h = (b-a)/4;
switch n
    case 1
out = (7*f(t,xval(1))*bval(5) + 32*f(t,xval(2))*bval(4) +...
              12*f(t,xval(3))*bval(3) + ...
         32*f(t,xval(4))*bval(2) + 7*f(t,xval(5))*bval(1))*(2*h/45);
    case 2
        out = (7*f(t,xval(1))*bval(1) + 32*f(t,xval(2))*bval(2) +...
              12*f(t,xval(3))*bval(3) + ...
         32*f(t,xval(4))*bval(4) + 7*f(t,xval(5))*bval(5))*(2*h/45);
end
     
% Booles rule has h^7 order error.