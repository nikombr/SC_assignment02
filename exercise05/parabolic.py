
import numpy as torch
import torch 
import torch.nn as nn
from utils import *

class Parabolic:

    def boundary_condition(self, x: torch.Tensor, t: torch.Tensor, epsilon: float = 0.1) -> torch.Tensor: # gL and gR
        self.epsilon = epsilon
        return self.u_fun(x,t)

    def initial_condition(self,x: torch.Tensor, t: torch.Tensor) -> torch.Tensor: # eta
        return self.u_fun(x,t)
    
    def u_fun(self,x,t): #this is exactly as in matlab, but it seems the expression doesn't match what is stated in overleaf for boundary-conditions
        alpha = [1, 4, 16]
        N = 3
        a = 1
        b = 0
        epsilon = 0.1
        u = 0.0
        for n in range(N):
            u += torch.exp(-epsilon * alpha[n]**2 * t) * (a*torch.cos(alpha[n]*x) + b*torch.sin(alpha[n]*x))
        return u
    
    def pde(self,model,x,t):
        return dfdt(model,x,t,order=1) - self.epsilon*dfdx(model,x,t,order=2)
