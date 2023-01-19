
# Numpy is used to preform vectorized operations on the solution arrays
import numpy as np

# Uses the finite difference method in 2 dimensions to solve the linear convection equation
# on a square array.
# Note: The number of iterative steps in x-direction is equal to the number of steps in the
# y-direction and taking the argument `sargs`. The time steps takes the argument `targs` with
# step size denoted by `dt`.
# Note: dx = dy. Denoted by the parameter `dx`.
# Note: `convection` is the convection coefficent or velocity field and is constant or uniform.
# Note: `init` is the initial condition.
class convection_2dims:
    '''
    ## Convection equation in 2 dimensions
    
    ### PARAMETERS:

    1. `targs: int` number of time steps.
    2. `sargs: int` number of spacial steps and is the same for both x-direction and y-direction.
    3. `dt: float` distance between time steps.
    4. `dx: float` distance between spacial steps.
    5. `convection: float` diffusion coefficent.
    6. `init: np.ndarray` initial condition.

    ### Methods:

    1. Solve the convection equation:

    >>>`convection_2dims.solve(self, boundary: list)`

    '''
    def __init__(self, targs: int, sargs: int, dt: float, dx: float, convection: float, init: np.ndarray) -> None:
        self.targs = targs
        self.sargs = sargs
        self.dt = dt
        self.dx = dx
        self.convection = convection
        self.init = init

    # Solves the 2d linear convection equation
    def solve(self) -> np.ndarray:
        '''
        ## Solves the Convection equation in 2 dimensions:

        Uses the finite differentce method.
        '''
        temp = []
        temp.append(self.init)
        
        for iargs in range(0, self.targs):

            # sets up next arrage of zeros in solution list
            next_array = np.zeros((self.sargs, self.sargs))
            temp.append(next_array)

            # Finite difference method in 2 dimensions
            temp[iargs + 1] = temp[iargs] + ((self.convection*self.dt)/(2*self.dx))*(
            np.roll(temp[iargs], 1, axis=1) - np.roll(temp[iargs], -1, axis=1) + 
            np.roll(temp[iargs], 1, axis=0) - np.roll(temp[iargs], -1, axis=0)
            )
            
        return temp

# Uses the finite difference method in 1 dimensions to solve the linear convection equation.
# Note: The number of iterative steps in x-direction takes the argument `sargs`. 
# The time steps takes the argument `targs` with step size denoted by `dt`.
# Spacial step size denoted by the parameter `dx`.
# Note: `convection` is the convection coefficent or velocity field and is constant or uniform.
# Note: `init` is the initial condition.
class convection_1dims:
    '''
    ## Convection equation in 1 dimensions
    
    ### PARAMETERS:

    1. `targs: int` number of time steps.
    2. `sargs: int` number of spacial steps.
    3. `dt: float` distance between time steps.
    4. `dx: float` distance between spacial steps.
    5. `convection: float` diffusion coefficent.
    6. `init: np.ndarray` initial condition.

    ### Methods:

    1. Solve the convection equation:

    >>>`convection_2dims.solve(self, boundary: list)`

    '''
    def __init__(self, targs: int, sargs: int, dt: float, dx: float, convection: float, init: np.ndarray) -> None:
        self.targs = targs
        self.sargs = sargs
        self.dt = dt
        self.dx = dx
        self.convection = convection
        self.init = init

    # Solves the 1d linear convection equation
    def solve(self) -> np.ndarray:
        '''
        ## Solves the Convection equation in 1 dimensions:

        Uses the finite differentce method.
        '''
        temp = []
        temp.append(self.init)
        
        for iargs in range(0, self.targs):

            # sets up next arrage of zeros in solution list
            next_array = np.zeros(self.sargs)
            temp.append(next_array)

            # Finite difference method in 2 dimensions
            temp[iargs + 1] = temp[iargs] + ((self.Convection*self.Deltat)/(2*self.Deltax))*(
            np.roll(temp[iargs], 1) - np.roll(temp[iargs], -1)
            )
            
        return temp

