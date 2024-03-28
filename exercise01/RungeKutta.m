function [T, Y, H, E, numStep] = RungeKutta(fun,y0,t0,tend,h0,reps,aeps,varargin)
% reps: relative error tolerance
% aeps: absolute error tolerance
% t0:   start time
% tend: end time
% y0:   initial value
% fun:  function handle
% Initialize
t = t0;
y = y0;
h = h0;
% Save data
T = [t];
Y = [y];
H = [h];
E = [];
numStep = 0;
facmin = 0.1; % Maximum decrease factor
facmax = 2.0;%5.0; % Maximum increase factor
% Iterate through time
while t < tend
    % If less than h away from end time - set h to be the time difference
    if (t + h > tend)
        h = tend - t;
    end
    % Do step
    [t, y, error] = RungeKuttaStep(fun,t,y,h,varargin{:});
    % Compute tol for step update
    tol = reps * norm(y) + aeps;
    % Asymptotic step size controller
    change = (tol/error)^(1/3); 
    %h = change*h;
    h = max(facmin, min(change,facmax))*h;
    % Save data
    T = [T t];
    Y = [Y y];
    H = [H h];
    E = [E error];
    % Count number of steps
    numStep = numStep + 1;
end