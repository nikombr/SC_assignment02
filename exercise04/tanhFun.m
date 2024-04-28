function u = tanhFun(x,t,epsilon)
u = -tanh((x + 1/2 - t)/(2*epsilon)) + 1;