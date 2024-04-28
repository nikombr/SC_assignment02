
import numpy as torch
import torch 
import torch.nn as nn
from utils import *

class Parabolic:

    def boundary_condition(self, x: torch.Tensor, t: torch.Tensor, epsilon: float = 0.1) -> torch.Tensor: # gL and gR
        if x[0] == -1:
            print("Found left boundary!")
            #g = torch.exp(-epsilon*t)*torch.cos(-1) + 4*torch.exp(-16*epsilon*t)*torch.cos(-4) + 16*torch.exp(-256*epsilon*t)*torch.cos(-16)
        elif x[0] == 1: 
            print("Found right boundary!")
            #g = torch.exp(-epsilon*t)*torch.cos(1) + 4*torch.exp(-16*epsilon*t)*torch.cos(4) + 16*torch.exp(-256*epsilon*t)*torch.cos(16)
        else:
            raise NotImplementedError(f"boundary_conditions received wrong boundary grid... x[0] = {x[0]}")
        self.epsilon = epsilon
        g = torch.exp(-epsilon*t)*torch.cos(1*x) + 4*torch.exp(-16*epsilon*t)*torch.cos(4*x) + 16*torch.exp(-256*epsilon*t)*torch.cos(16*x)
        return g 

    def initial_condition(self,x: torch.Tensor, t: torch.Tensor) -> torch.Tensor: # eta
        assert t[0] == 0.0, f"Called eta boundary condition for wrong t value: {t} (should be 0.0)"
        eta = torch.cos(x)+4*torch.cos(4*x) + 16*torch.cos(16*x)
        return eta
    
    def u_fun(x,t): #this is exactly as in matlab, but it seems the expression doesn't match what is stated in overleaf for boundary-conditions
        alpha = [1, 4, 16]
        N = 3
        a = 1
        b = 0
        epsilon = 0.1
        u = 0.0
        for n in range(N):
            #first = torch.exp(-epsilon*alpha[n]**2*T)
            #second = a*torch.cos(alpha[n]*X)
            #third = b*torch.sin(alpha[n]*X)
            #print(first.shape, second.shape, third.shape)
            u += torch.exp(-epsilon*alpha[n]**2*t)*(a*torch.cos(alpha[n]*x) + b*torch.sin(alpha[n]*x))
        return u
    
    def pde(self,model,x,t):
        return dfdt(model,x,t,order=1) - self.epsilon*dfdx(model,x,t,order=2)
