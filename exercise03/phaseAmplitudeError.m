function [amplitude_error, phase_error] = phaseAmplitudeError(utrue,U,x,t)

amplitude_error = zeros(length(t)-1,1)

for i = 2:length(t)
    
    amplitude_error(i-1) = max(utrue(i,:)) - max(U(i,:))


end

figure;
plot(amplitude_error'./t(2:end))

amplitude_error = mean(amplitude_error'./t(2:end))

phase_error = 0