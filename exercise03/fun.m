function u = fun(x,t,a)

if nargin == 2
    a = 0.5;
end

u = sin(2*pi*(x-a*t));