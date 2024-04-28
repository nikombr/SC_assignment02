import numpy as np 
import torch 
import torch.nn as nn
import matplotlib.pyplot as plt 

class PINN(nn.Module):
    def __init__(self,hidden_dimensions: list = [2,128,1],activation_fn = nn.Tanh()):
        super(PINN,self).__init__()
        self.layers = nn.ModuleList([])
        self.activation_fn = activation_fn
        for i in range(len(hidden_dimensions)-1):
            self.layers.append(nn.Linear(hidden_dimensions[i],hidden_dimensions[i+1]))
            if(i < len(hidden_dimensions)-2): #no activation for final layer 
                self.layers.append(self.activation_fn)
        
        self.layers = nn.Sequential(*self.layers)

    def forward(self,x,t):
        #print_shape(x)
        #print_shape(t)
        x_stack = torch.cat([x,t],dim=1)
        return self.layers(x_stack)
    
#order is order of differentiation
def df(output: torch.Tensor, input_var: torch.Tensor, order: int = 1) -> torch.Tensor:
    """Compute neural network derivative with respect to input features using PyTorch autograd engine"""
    df_value = output      # <-- we directly take the output of the NN
    for _ in range(order):
        df_value = torch.autograd.grad(
            df_value,
            input_var,
            grad_outputs=torch.ones_like(input_var),
            create_graph=True,
            retain_graph=True,
        )[0]
    return df_value

def dfdt(model: PINN, x: torch.Tensor, t: torch.Tensor, order: int = 1):
    """Derivative with respect to the time variable of arbitrary order"""
    f_value = model(x, t)
    return df(f_value, t, order=order)

def dfdx(model: PINN, x: torch.Tensor, t: torch.Tensor, order: int = 1):
    """Derivative with respect to the spatial variable of arbitrary order"""
    f_value = model(x, t)
    return df(f_value, x, order=order)