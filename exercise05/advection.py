
import numpy as torch
import torch 
import torch.nn as nn
from utils import *

class Advection:

    def __init__(self, epsilon: float = 0.01/torch.pi) -> None:
        self.epsilon  = epsilon

    def boundary_condition(self, x: torch.Tensor, t: torch.Tensor, epsilon: float = 0.1) -> torch.Tensor: # gL and gR
        return torch.zeros_like(x)

    def initial_condition(self,x: torch.Tensor, t: torch.Tensor) -> torch.Tensor: # eta
        assert t[0] == 0.0, f"Called eta boundary condition for wrong t value: {t} (should be 0.0)"
        eta = -torch.sin(torch.pi*x)
        return eta
    
    def pde(self,model,x,t):
        return dfdt(model,x,t,order=1) + dfdx(model,x,t,order=0)*dfdx(model,x,t,order=1) - self.epsilon*dfdx(model,x,t,order=2)
