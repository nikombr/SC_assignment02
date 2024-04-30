
from utils import *
from parabolic import *
from hyperbolic import *
from advection import *
from torch.optim import Adam 
from tqdm import tqdm
import inspect
from scipy.io import savemat
from odesolver import *
import sys

runtype = int(sys.argv[1])

N = 1000
M = 1000

if runtype == 1: # Parabolic

    name = "parabolic"
    parabolic = Parabolic()
    tmax = 0.5

    # Test of parameters
    
    hn = [20, 40, 80]

    for i, h in enumerate(hn):
        H = [[2,h,1], [2,h,h,h,1], [2,h,h,h,h,h,1]]

        for j in range(3):

            torch.manual_seed(1)

            model = ODESolver(hidden_dimensions = H[j],
                            activation_fn = nn.Tanh(),
                            N = N,
                            M = M,
                            tmax = tmax,
                            verbose = False,
                            pde = parabolic)

            N_EPOCHS = 10000
            loss = model.train(N_EPOCHS)

            U = model.predict().cpu().detach().numpy().reshape([N,M]).T
            X = model.X.detach().cpu().numpy().reshape([N,M]).T
            T = model.T.detach().cpu().numpy().reshape([N,M]).T
            x = model.x.detach().cpu().numpy()
            t = model.t.detach().cpu().numpy()
            mdic = {"U": U, "T": T, "X": X, "loss": loss, "x": x, "t": t}

            savemat(f"results/{name}_{i}_{j}.mat", mdic)
    
    # Convergence test

    M = 10000

    H = [2,40,40,40,40,40,1];

    for N in np.round(np.linspace(10,200,10)):
        
        N = int(N)

        torch.manual_seed(1)

        model = ODESolver(hidden_dimensions = H,
                        activation_fn = nn.Tanh(),
                        N = N,
                        M = M,
                        tmax = tmax,
                        verbose = False,
                        pde = parabolic)

        N_EPOCHS = 10000
        loss = model.train(N_EPOCHS)

        U = model.predict().cpu().detach().numpy().reshape([N,M]).T
        X = model.X.detach().cpu().numpy().reshape([N,M]).T
        T = model.T.detach().cpu().numpy().reshape([N,M]).T
        x = model.x.detach().cpu().numpy()
        t = model.t.detach().cpu().numpy()
        mdic = {"U": U, "T": T, "X": X, "loss": loss, "x": x, "t": t}

        savemat(f"results/{name}_convergence_spatial_N_{N}.mat", mdic)

    for N in np.round(np.linspace(10,200,10)):
        N = int(N)
        h = 2/N
        k = h**2
        M = int(np.ceil(tmax/k))

        torch.manual_seed(1)

        model = ODESolver(hidden_dimensions = H,
                        activation_fn = nn.Tanh(),
                        N = N,
                        M = M,
                        tmax = tmax,
                        verbose = False,
                        pde = parabolic)

        N_EPOCHS = 10000
        loss = model.train(N_EPOCHS)

        U = model.predict().cpu().detach().numpy().reshape([N,M]).T
        X = model.X.detach().cpu().numpy().reshape([N,M]).T
        T = model.T.detach().cpu().numpy().reshape([N,M]).T
        x = model.x.detach().cpu().numpy()
        t = model.t.detach().cpu().numpy()
        mdic = {"U": U, "T": T, "X": X, "loss": loss, "x": x, "t": t}

        savemat(f"results/{name}_convergence_combined_N_{N}.mat", mdic)


elif runtype == 2: # Hyperbolic

    name = "hyperbolic"

    hyperbolic = Hyperbolic()

    H = [2,40,40,40,40,40,1];

    tmax = 2.0
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


    model = ODESolver(hidden_dimensions = H,
                    activation_fn = nn.Tanh(),
                    N = N,
                    M = M,
                    tmax = tmax,
                    verbose = False,
                    pde = hyperbolic)

    N_EPOCHS = 10000
    loss = model.train(N_EPOCHS)

    U = model.predict().cpu().detach().numpy().reshape([N,M]).T
    X = model.X.detach().cpu().numpy().reshape([N,M]).T
    T = model.T.detach().cpu().numpy().reshape([N,M]).T
    x = model.x.detach().cpu().numpy()
    t = model.t.detach().cpu().numpy()
    mdic = {"U": U, "T": T, "X": X, "loss": loss, "x": x, "t": t}

    # Predict outside of traning time
    # Create grid over the whole domain
    x_pred   = torch.linspace(-1, 1, steps = N, requires_grad=True).to(device)
    t_pred   = torch.linspace(tmax   , 2*tmax, steps = M, requires_grad=True).to(device)
    grid   = torch.meshgrid(x_pred, t_pred, indexing="ij")
    X_pred, T_pred = grid[0].flatten().reshape(-1,1).to(device), grid[1].flatten().reshape(-1,1).to(device)
    U_pred = model.predict(x=X_pred,t=T_pred).cpu().detach().numpy().reshape([N,M]).T

    mdic = {"U": U, "T": T, "X": X, "loss": loss, "x": x, "t": t, "U_pred": U_pred, "t_pred": t_pred.detach().cpu().numpy()}

elif runtype == 3: # Advection

    name = "advection"

    advection = Advection()

    tmax = 1.6037/torch.pi # Time where we want to estiamte the derivative

    H = [2,20,20,20,20,20,1];
    H = [2,40,40,40,40,40,1];

    model = ODESolver(hidden_dimensions = H,
                    activation_fn = nn.Tanh(),
                    N = N,
                    M = M,
                    tmax = tmax,
                    verbose = False,
                    pde = advection)

    N_EPOCHS = 10000
    loss = model.train(N_EPOCHS)

    U = model.predict().cpu().detach().numpy().reshape([N,M]).T
    X = model.X.detach().cpu().numpy().reshape([N,M]).T
    T = model.T.detach().cpu().numpy().reshape([N,M]).T
    x = model.x.detach().cpu().numpy()
    t = model.t.detach().cpu().numpy()
    mdic = {"U": U, "T": T, "X": X, "loss": loss, "x": x, "t": t}


    savemat(f"results/{name}.mat", mdic)