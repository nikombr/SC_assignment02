
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

    model = ODESolver(hidden_dimensions = [2,20,20,20,20,20,1],
                    activation_fn = nn.Tanh(),
                    N = N,
                    M = M,
                    tmax = 1.0,
                    verbose = False,
                    pde = parabolic)

    N_EPOCHS = 3000
    loss = model.train(N_EPOCHS)

elif runtype == 2: # Hyperbolic

    name = "hyperbolic"

    hyperbolic = Hyperbolic()

    model = ODESolver(hidden_dimensions = [2,20,20,20,20,20,1],
                    activation_fn = nn.Tanh(),
                    N = N,
                    M = M,
                    tmax = 2.0,
                    verbose = False,
                    pde = hyperbolic)

    N_EPOCHS = 3000
    loss = model.train(N_EPOCHS)

elif runtype == 3: # Advection

    name = "advection"

    advection = Advection()

    tmax = 1.6037/torch.pi # Time where we want to estiamte the derivative

    model = ODESolver(hidden_dimensions = [2,20,20,20,20,20,1],
                    activation_fn = nn.Tanh(),
                    N = N,
                    M = M,
                    tmax = tmax,
                    verbose = False,
                    pde = advection)

    N_EPOCHS = 3000
    loss = model.train(N_EPOCHS)


U = model.predict().cpu().detach().numpy().reshape([N,M]).T
X = model.X.detach().cpu().numpy().reshape([N,M]).T
T = model.T.detach().cpu().numpy().reshape([N,M]).T
mdic = {"U": U, "T": T, "X": X, "loss": "loss"}
savemat(f"results/{name}.mat", mdic)