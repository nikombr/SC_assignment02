function [time, funEval, stepSize, error, numStep] = RungeKutta(fun,y0,t0,tend,h0,reps,aeps,varargin)
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
time = [t]; funEval = [y]; stepSize = [h]; error = [];
numStep = 0;
%facmin = 0.1; % Maximum decrease factor
%facmax = 2.0;%5.0; % Maximum increase factor
% Iterate through time
while t < tend
    % If less than h away from end time - set h to be the time difference
    if (t + h > tend)
        h = tend - t;
    end
    % Do step
    [t, y, e] = RungeKuttaStep(fun,t,y,h,varargin{:});
    % Compute tol for step update
    tol = reps * norm(y) + aeps;
    % Asymptotic step size controller
    h = (tol/e)^(1/3)*h;
    % Save data
    time = [time t]; funEval = [funEval y]; stepSize = [stepSize h]; error = [error e];
    % Count number of steps
    numStep = numStep + 1;
end