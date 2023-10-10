
# Numpy is used to preform vectorized operations on the solution arrays
import numpy as np

# Diffusion Convection model in 2 dimensions
class convection_diffusion_2dims:
    '''
    ## Diffusion equation in 2 dimensions
    
    ### PARAMETERS:

    1. `targs: int` number of time steps.
    2. `sargs: int` number of spacial steps and is the same for both x-direction and y-direction.
    3. `dt: float` distance between time steps.
    4. `dx: float` distance between spacial steps.
    5. `diffusion: float` diffusion coefficient.
    6. `convection: float` convection coeffient.
    7. `init: np.ndarray` initial condition.

    ### Methods:

    1. Solve the diffusion equation under Dirichlet Boundary Conditions:

    >>>`diffusion_2dims.solve_Dirichlet(self, boundary: list)`

    2. Solve the diffusion equation under Neumann Boundary Conditions:

    >>>`diffusion_2dims.solve_Neumann(self, boundary: list)` 


    '''

    def __init__(self, targs: int, sargs: int, dt: float, dx: float, diffusion: float, convection: float, init: np.ndarray):
        self.targs = targs
        self.sargs = sargs
        self.dt = dt
        self.dx = dx
        self.diffusion = diffusion
        self.convection = convection
        self.init = init

    def solve_Dirichlet(self, boundary: list):
        """
        Solves the 2 dimensional Diffusion Convection Equation under Dirichlet Boundary Conditions.
        """

        if len(boundary) != 4:
            raise ValueError(f'Boundary list should have a lenght of 4 but boundary list is given with lenght {len(boundary)}')  
                
        solution = []
        solution.append(self.init)
        for iargs in range(0, self.targs):

            # sets up next arrage of zeros in solution list
            next_array = np.zeros((self.sargs, self.sargs))
            solution.append(next_array)

            # Finite difference method in 2 dimensions
            solution[iargs + 1] = solution[iargs] + ((self.diffusion*self.dt)/self.dx**2)*(
            np.roll(solution[iargs], 1, axis=1) + np.roll(solution[iargs], -1, axis=1) + 
            np.roll(solution[iargs], 1, axis=0) + np.roll(solution[iargs], -1, axis=0) - 4*solution[iargs]
            ) + ((self.convection*self.dt)/(2*self.dx))*(
            np.roll(solution[iargs], 1, axis=1) - np.roll(solution[iargs], -1, axis=1) + 
            np.roll(solution[iargs], 1, axis=0) - np.roll(solution[iargs], -1, axis=0)
            )

            # Dirichlet boundary conditions
            solution[iargs][:, 0] = boundary[0]
            solution[iargs][:, -1] = boundary[1]
            solution[iargs][0, :] = boundary[2]
            solution[iargs][-1, :] = boundary[3]
                
        return _pde('2 Dimensions', 'Dirichlet', self.diffusion, self.convection, np.array(solution), boundary, self.init, self.dt, self.dx, [self.sargs, self.targs])
            

    def solve_Neumann(self, boundary_flux: list):
        """
        Solves the 2 dimensional Diffusion Convection Equation under Neumann Boundary Conditions.
        """

        if len(boundary_flux) != 4:
            raise ValueError(f'Boundary flux list should have a lenght of 4 but boundary_flux list is given with lenght {len(boundary_flux)}')  
                
        solution = []
        solution.append(self.init)
        for iargs in range(0, self.targs):

            # sets up next arrage of zeros in solution list
            next_array = np.zeros((self.sargs, self.sargs))
            solution.append(next_array)

            # Finite difference method in 2 dimensions
            solution[iargs + 1] = solution[iargs] + ((self.diffusion*self.dt)/self.dx**2)*(
            np.roll(solution[iargs], 1, axis=1) + np.roll(solution[iargs], -1, axis=1) + 
            np.roll(solution[iargs], 1, axis=0) + np.roll(solution[iargs], -1, axis=0) - 4*solution[iargs]
            ) + ((self.convection*self.dt)/(2*self.dx))*(
            np.roll(solution[iargs], 1, axis=1) - np.roll(solution[iargs], -1, axis=1) + 
            np.roll(solution[iargs], 1, axis=0) - np.roll(solution[iargs], -1, axis=0)
            )

            # Neumann boundary conditions
            solution[iargs][:, 0] = solution[iargs][:, 1] - boundary_flux[0]*self.dx
            solution[iargs][:, -1] = solution[iargs][:, -2] + boundary_flux[1]*self.dx
            solution[iargs][0, :] = solution[iargs][1, :] - boundary_flux[2]*self.dx
            solution[iargs][-1, :] = solution[iargs][-2, :] + boundary_flux[3]*self.dx
                
        return _pde('2 Dimensions', 'Neumann', self.diffusion, self.convection, np.array(solution), boundary_flux, self.init, self.dt, self.dx, [self.sargs, self.targs])
            
# Diffusion Convection model in 1 dimension.
class convection_diffusion_1dims:
    '''
    ## Diffusion equation in 1 dimension
    
    ### PARAMETERS:

    1. `targs: int` number of time steps.
    2. `sargs: int` number of spacial steps and is the same for both x-direction and y-direction.
    3. `dt: float` distance between time steps.
    4. `dx: float` distance between spacial steps.
    5. `diffusion: float` diffusion coefficient.
    6. `convection: float` convection coeffient.
    7. `init: np.ndarray` initial condition.

    ### Methods:

    1. Solve the diffusion equation under Dirichlet Boundary Conditions:

    >>>`diffusion_2dims.solve_Dirichlet(self, boundary: list)`

    2. Solve the diffusion equation under Neumann Boundary Conditions:

    >>>`diffusion_2dims.solve_Neumann(self, boundary: list)` 


    '''

    def __init__(self, targs: int, sargs: int, dt: float, dx: float, diffusion: float, convection: float, init: np.ndarray):
        self.targs = targs
        self.sargs = sargs
        self.dt = dt
        self.dx = dx
        self.diffusion = diffusion
        self.convection = convection
        self.init = init

    def solve_Dirichlet(self, boundary: list):
        """
        Solves the 1 dimensional Diffusion Convection Equation under Dirichlet Boundary Conditions.
        """

        if len(boundary) != 2:
            raise ValueError(f'Boundary list should have a lenght of 2 but boundary list is given with lenght {len(boundary)}')  
        
        solution = []
        solution.append(self.init)
        for iargs in range(0, self.targs):

            # sets up next arrage of zeros in solution list 
            next_array = np.zeros(self.sargs)
            solution.append(next_array)

            # Finite difference method in 1 dimension
            solution[iargs + 1] = solution[iargs] + ((self.diffusion*self.dt)/(self.dx**2))*(
            np.roll(solution[iargs], 1) + np.roll(solution[iargs], -1) - 2*solution[iargs]
            ) + ((self.convection*self.dt)/(2*self.dx))*(
            np.roll(solution[iargs], 1) - np.roll(solution[iargs], -1)
            )

            # Dirichlet boundary conditions
            solution[iargs][ 0] = boundary[0]
            solution[iargs][-1] = boundary[1]
                
        return _pde('1 Dimension', 'Dirichlet', self.diffusion, self.convection, np.array(solution), boundary, self.init, self.dt, self.dx, [self.sargs, self.targs])


    def solve_Neumann(self, boundary_flux: list):
        """
        Solves the 1 dimensional Diffusion Convection Equation under Neumann Boundary Conditions.
        """

        if len(boundary_flux) != 2: 
            raise ValueError(f'Boundary flux list should have a lenght of 2 but boundary_flux list is given with lenght {len(boundary_flux)}') 
                
        solution = []
        solution.append(self.init)
        for iargs in range(0, self.targs):

            # sets up next arrage of zeros in solution list 
            next_array = np.zeros(self.sargs)
            solution.append(next_array)

            # Finite difference method in 1 dimension
            solution[iargs + 1] = solution[iargs] + ((self.diffusion*self.dt)/(self.dx**2))*(
            np.roll(solution[iargs], 1) + np.roll(solution[iargs], -1) - 2*solution[iargs]
            ) + ((self.convection*self.dt)/(2*self.dx))*(
            np.roll(solution[iargs], 1) - np.roll(solution[iargs], -1)
            )

            # Dirichlet boundary conditions
            solution[iargs][0] = solution[iargs][1] - boundary_flux[0]*self.dx
            solution[iargs][-1] = solution[iargs][-2] + boundary_flux[1]+self.dx
                
        return _pde('1 Dimension', 'Neumann', self.diffusion, self.convection, np.array(solution), boundary_flux, self.init, self.dt, self.dx, [self.sargs, self.targs])



class _pde:
    def __init__(self, pde_type: str, pde_boundary_type: str, diffusion: float, convection: float, solution: np.ndarray, boundary_conditions: list, initial_state: np.ndarray, dt: float, dx: float, step_size: list):
        self.pde_type = pde_type
        self.pde_boundary_type = pde_boundary_type
        self.diffusion = diffusion
        self.convection = convection
        self.solution = solution
        self.boundary_conditions = boundary_conditions
        self.initial_state = initial_state
        self.dt = dt
        self.dx = dx
        self.step_size = step_size

    def __repr__(self):
        return f'Diffusion in {self.pde_type} with {self.pde_boundary_type} Conditions\nDiffusion Coefficient: {self.diffusion}\nConvection Coefficient: {self.convection}\nSolution: {self.solution}\nBoundary Conditions {self.boundary_conditions}\nInitial State: {self.initial_state}\ndt: {self.dt}\ndx = dy: {self.dx}\nNumber of steps in x and y: {self.step_size[0]}\nNumber of time steps {self.step_size[1]}'               