
import numpy as torch
import torch 
import torch.nn as nn
from utils import *

class Hyperbolic:

    def __init__(self, a: float = 0.5) -> None:
        self.a = a
        self.boundary_condition = False

    def initial_condition(self,x: torch.Tensor, t: torch.Tensor) -> torch.Tensor: # eta
        eta = torch.sin(2*torch.pi*x)
        return eta
    
    def u_fun(self,x,t): #this is exactly as in matlab, but it seems the expression doesn't match what is stated in overleaf for boundary-conditions
        u = torch.sin(2*torch.pi*(x-self.a*t))
        return u
    
    def pde(self,model,x,t):
        return dfdt(model,x,t,order=1) + self.a*dfdx(model,x,t,order=1)