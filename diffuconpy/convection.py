
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
    def __init__(self, targs: int, sargs: int, dt: float, dx: float, convection: float, init: np.ndarray):
        self.targs = targs
        self.sargs = sargs
        self.dt = dt
        self.dx = dx
        self.convection = convection
        self.init = init

    # Solves the 2d linear convection equation
    def solve(self):
        '''
        ## Solves the Convection equation in 2 dimensions:

        Uses the finite difference method.
        '''
        solution = []
        solution.append(self.init)
        
        for iargs in range(0, self.targs):

            # sets up next arrage of zeros in solution list
            next_array = np.zeros((self.sargs, self.sargs))
            solution.append(next_array)

            # Finite difference method in 2 dimensions
            solution[iargs + 1] = solution[iargs] + ((self.convection*self.dt)/(2*self.dx))*(
            np.roll(solution[iargs], 1, axis=1) - np.roll(solution[iargs], -1, axis=1) + 
            np.roll(solution[iargs], 1, axis=0) - np.roll(solution[iargs], -1, axis=0)
            )
            
        return _pde('2 Dimensions', np.array(solution), self.init, self.dt, self.dx, [self.sargs, self.targs])

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

    >>>`convection_1dims.solve(self, boundary: list)`

    '''
    def __init__(self, targs: int, sargs: int, dt: float, dx: float, convection: float, init: np.ndarray):
        self.targs = targs
        self.sargs = sargs
        self.dt = dt
        self.dx = dx
        self.convection = convection
        self.init = init

    # Solves the 1d linear convection equation
    def solve(self):
        '''
        ## Solves the Convection equation in 1 dimensions:

        Uses the finite difference method.
        '''
        solution = []
        solution.append(self.init)
        
        for iargs in range(0, self.targs):

            # sets up next arrage of zeros in solution list
            next_array = np.zeros(self.sargs)
            solution.append(next_array)

            # Periodic Boundary Conditions
            solution[iargs + 1][0] = solution[iargs][0] +  ((self.convection*self.dt)/(2*self.dx))*(
                solution[iargs][1] - solution[iargs][self.sargs-1]
            )

            solution[iargs + 1][self.sargs-1] = solution[iargs][self.sargs-1] + ((self.convection*self.dt)/(2*self.dx))*(
                solution[iargs][0] - solution[iargs][self.sargs-2]
            )

            #solution[iargs + 1][0] = solution[iargs + 1][-1]

            # Finite difference method in 2 dimensions
            solution[iargs + 1] = solution[iargs] + ((self.convection*self.dt)/(2*self.dx))*(
            np.roll(solution[iargs], 1) - np.roll(solution[iargs], -1)
            )
            
        return _pde('1 Dimension', np.array(solution), self.init, self.dt, self.dx, [self.sargs, self.targs])


class _pde:
    def __init__(self, pde_type: str, solution: np.ndarray, initial_state: np.ndarray, dt: float, dx: float, step_size: list):
        self.pde_type = pde_type
        self.solution = solution
        self.initial_state = initial_state
        self.dt = dt
        self.dx = dx
        self.step_size = step_size

    def __repr__(self):
        return f'Convection in {self.pde_type}\nSolution: {self.solution}\nInitial State: {self.initial_state}\ndt: {self.dt}\ndx = dy: {self.dx}\nNumber of steps in x and y: {self.step_size[0]}\nNumber of time steps {self.step_size[1]}'
        