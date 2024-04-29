
from utils import *
from torch.optim import Adam 
from tqdm import tqdm
import torch.nn as nn

class ODESolver:
    def __init__(self,
                 hidden_dimensions: list = None,
                 activation_fn = nn.Tanh(),
                 tmax: float  = 1, # Time is in [0, tmax]
                 M: int = 100, # Number of time points
                 N: int = 100, # Number of spatial points
                 xmin: float = -1, 
                 xmax: float  = 1, 
                 verbose = False,
                 pde = False,
                 learning_rate: float = 1e-2
                ):

        # Check for cuda

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        self.verbose = verbose
        self.model = PINN(hidden_dimensions = hidden_dimensions,activation_fn = activation_fn)
        
        # Create grid over the whole domain
        x   = torch.linspace(xmin, xmax, steps = N, requires_grad=True).to(device)
        t   = torch.linspace(0   , tmax, steps = M, requires_grad=True).to(device)
        grid   = torch.meshgrid(x, t, indexing="ij")

        # Save grid to class
        self.x, self.t = x.reshape(-1,1).to(device), t.reshape(-1,1).to(device)
        self.X, self.T = grid[0].flatten().reshape(-1,1).to(device), grid[1].flatten().reshape(-1,1).to(device)

        # Get pde information
        self.initial_condition  = pde.initial_condition
        self.boundary_condition = pde.boundary_condition
        self.pde                = pde.pde

        # Set place of initial and boundary conditions
        self.t_initial = (torch.ones_like(x, requires_grad = True)*t[0]).reshape(-1,1).to(device)
        self.x_boundary_gL = (torch.ones_like(t, requires_grad = True)*x[0]).reshape(-1,1).to(device)
        self.x_boundary_gR = (torch.ones_like(t, requires_grad = True)*x[-1]).reshape(-1,1).to(device)

        # Get the value in initial and boundary conditions
        self.u_initial     = self.initial_condition(x = self.x, t = self.t_initial).to(device)
        if self.boundary_condition is not False:
            self.u_boundary_gL = self.boundary_condition(x = self.x_boundary_gL, t = self.t).to(device)
            self.u_boundary_gR = self.boundary_condition(x = self.x_boundary_gR, t = self.t).to(device)
                
        #instantiate loss criterion: 
        self.criterion = torch.nn.MSELoss()
        
        #setup optimization
        self.optimizer = Adam(self.model.parameters(), lr=learning_rate, betas=(0.9,0.999), weight_decay=0, amsgrad=False)
        

    def predict_compute_losses(self):

        # Predict with model
        u_pred_initial      = self.model(self.x,self.t_initial)
        u_pred_boundary_gL  = self.model(self.x_boundary_gL,self.t)
        u_pred_boundary_gR  = self.model(self.x_boundary_gR,self.t)
        u_pred_pde          = self.pde(self.model,self.X,self.T)

        # Loss from initial condition
        initial_loss        = self.criterion(u_pred_initial,self.u_initial)

        # Loss from boundary condition
        if self.boundary_condition is False: # We need the two boundaries to be equal
            boundary_loss       = self.criterion(u_pred_boundary_gL, u_pred_boundary_gR)
        else:
            boundary_gL_loss    = self.criterion(u_pred_boundary_gL,self.u_boundary_gL)
            boundary_gR_loss    = self.criterion(u_pred_boundary_gR,self.u_boundary_gR)
            boundary_loss       = boundary_gL_loss + boundary_gR_loss

        # Loss from PDE
        pde_loss            = self.criterion(u_pred_pde,torch.zeros_like(u_pred_pde))

        # Total loss
        loss = initial_loss + boundary_loss + pde_loss
        

        return loss

    
    def train(self,nepochs=5000):
        self.model.train()
        loss_holder = []
        pbar = tqdm(range(nepochs))
        for i in pbar:
            loss = self.predict_compute_losses()
            self.optimizer.zero_grad()
            if(i % 5 == 0):
                pbar.set_description_str(f"epoch {i} | loss: {loss.item():.4f}")
            loss.backward(retain_graph=True)
            nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
            self.optimizer.step()
            loss_holder.append(loss.item())
        return loss_holder
    
    
    def predict(self,x=None,t=None):
        self.model.eval()
        #X = torch.stack(torch.meshgrid(x,t)).reshape(2,-1).T
        if (x == None or t == None):
            y_pred = self.model(x=self.X,t=self.T)
        else:
            y_pred = self.model(x,t)
        #y_pred = y_pred.reshape(len(x),len(t)).detach().cpu().numpy()
        self.model.train()
        return y_pred