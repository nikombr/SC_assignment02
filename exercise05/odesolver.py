
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
                 pde = False
                ):

        # Check for cuda

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        """
        arguments: self explanatory
            differentials: list of booleans encoding truth value for needed dx, dt, dxx, dtt.
            differential_equation: a lambda function calculating the network output and ground truth argument for the loss-criterion
        """
        #device = torch.device("cuda") if torch.cuda_is_available() else torch.device("cpu")
        self.verbose = verbose
        self.model = PINN(hidden_dimensions = hidden_dimensions,activation_fn = activation_fn)#PINN(hidden_dimensions=hidden_dimensions,activation_fn=activation_fn)
        
        # Create grid over the whole domain
        x   = torch.linspace(xmin, xmax, steps = N, requires_grad=True).to(device)
        t   = torch.linspace(0   , tmax, steps = M, requires_grad=True).to(device)
        grid   = torch.meshgrid(x, t, indexing="ij")

        # Save grid to class
        self.x, self.t = x.reshape(-1,1).to(device), t.reshape(-1,1).to(device)
        self.X, self.T = grid[0].flatten().reshape(-1,1).to(device), grid[1].flatten().reshape(-1,1).to(device)

        #self.X = torch.stack(torch.meshgrid(x, t)).reshape(2, -1).T

        # Get pde information
        self.initial_condition  = pde.initial_condition
        self.boundary_condition = pde.boundary_condition
        self.pde                = pde.pde

        #  
        #grid_gL   = torch.stack(torch.meshgrid(x[0],  t)).reshape(2, -1).T
        #grid_gR   = torch.stack(torch.meshgrid(x[-1], t)).reshape(2, -1).T
        #grid_eta  = torch.stack(torch.meshgrid(x,  t[0])).reshape(2, -1).T
        #self.X_train = torch.cat([grid_gL, grid_gR, grid_eta])
        #gL       = self.boundary_condition(grid_gL)
        #gR       = self.boundary_condition(grid_gR)
        #eta      = self.initial_condition(grid_eta)
        #self.y_train = torch.cat([gL, gR, eta]).unsqueeze(1)

        # Set place of initial and boundary conditions
        self.t_initial = (torch.ones_like(x, requires_grad = True)*t[0]).reshape(-1,1).to(device)
        self.x_boundary_gL = (torch.ones_like(t, requires_grad = True)*x[0]).reshape(-1,1).to(device)
        self.x_boundary_gR = (torch.ones_like(t, requires_grad = True)*x[-1]).reshape(-1,1).to(device)

        # Get the value in initial and boundary conditions
        self.u_initial     = self.initial_condition(x = self.x, t = self.t_initial).to(device)
        if self.boundary_condition is not False:
            self.u_boundary_gL = self.boundary_condition(x = self.x_boundary_gL, t = self.t).to(device)
            self.u_boundary_gR = self.boundary_condition(x = self.x_boundary_gR, t = self.t).to(device)

        #save everything we need to class
        #self.x, self.t = grids[0].flatten().reshape(-1,1), grids[1].flatten().reshape(-1,1)
        #self.t_initial = torch.ones_like(x_idx, requires_grad=True)*t_idx[0]
        #self.t_initial = self.t_initial.reshape(-1,1)
        #self.y_initial = initial_condition(self.x_idx) #should be made a function handle 
        #self.boundary_x0 = torch.ones_like(self.t_idx,requires_grad=True)*self.x[0]
        #self.boundary_x1 = torch.ones_like(self.t_idx,requires_grad=True)*self.x[-1]
        
        #todo: 
        #add boundary_y0 and boundary_y1 (currently working with initial condition f(x0,t)=f(x1,t)=0, but it should have rhs)
        #also note that we in this case have added an initial condition on f, but ALSO on dfdt(t=0)
                
        #instantiate loss criterion: 
        self.criterion = torch.nn.MSELoss()
        
        #setup optimization
        self.optimizer = Adam(self.model.parameters(), lr=1e-2, betas=(0.9,0.999), weight_decay=0.0, amsgrad=False)
        

    def predict_compute_losses(self):
        #pred_x0 = self.model(x = self.boundary_x0, t = self.t_idx) #predict at lower boundary
        #pred_x1 = self.model(x = self.boundary_x1, t = self.t_idx)

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
        
        #the way I'd rather want to do it later
        #initial_loss = self.criterion(pred_initial - self.y_initial)
        #initial_loss_df = self.criterion(dfdt(self.model,self.x_idx,self.t_initial,order=1))
        #pde_loss = self.criterion(df)
        
        #this is the magic - but realize it is almost exactly what you did earlier. Recheck your loss, recheck your dimensions. Your implementation seems fine
        #initial_loss_f = self.criterion(self.model(self.x_idx,self.t_initial),self.y_initial) #get mse-loss f(t=0( between initial condition anx y predicted at initial  
        #dfdt_ = dfdt(self.model,self.x_idx,self.t_initial,order=1) #needed for "pde loss"
        #initial_loss_df = self.criterion(dfdt_,torch.zeros_like(dfdt_))  #condition says it should be zero, so form a loss based on that
        #C=1.0 #constant for specific pde
        #pde_loss = self.criterion(dfdx(self.model,self.x,self.t,order=2),(1/C**2)*dfdt(self.model,self.x,self.t,order=2)) #pde-loss is basically left hand side minus right hand side, so either set target zero, or set output = lhs, target = rhs
        #loss = pde_loss + initial_loss_df + initial_loss_f + self.criterion(pred_x0,torch.zeros_like(pred_x0)) + self.criterion(pred_x1,torch.zeros_like(pred_x1)) #boundary conditions in this case are zero :)
        
        #this stuff also worked  if it is nicer to read:
        #initial_loss_f = self.model(self.x_idx,self.t_initial) - self.y_initial
        #initial_loss_df = dfdt(self.model,self.x_idx,self.t_initial,order=1)
        #C = 1.0
        #pde_loss = dfdx(self.model,self.x,self.t,order=2) - (1/C**2) * dfdt(self.model,self.x,self.t,order=2)
        #loss = \
        #pde_loss.pow(2).mean() + \
        #pred_x0.pow(2).mean() + \
        #pred_x1.pow(2).mean() + \
        #initial_loss_f.pow(2).mean() + \
        #initial_loss_df.pow(2).mean()

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