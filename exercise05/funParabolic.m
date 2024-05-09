function u = funParabolic(x,t)

alpha = [1 4 16];
N = 3;
a = 1;
b = 0;
epsilon = 0.1;

u = 0;

for n = 1:N
    u = u + exp(-epsilon*alpha(n)^2*t).*(a*cos(alpha(n)*x) + b*sin(alpha(n)*x));
end