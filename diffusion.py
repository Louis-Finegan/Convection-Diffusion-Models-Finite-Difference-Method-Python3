
# Numpy is used
import numpy as np

# Uses the finite difference method in 2 dimensions to solve the linear diffusion equation
# on a square array. Dirichlet boundary conditions are automaticly applied in the `solve_Dirichlet` method.
# Neumann boundary conditions can be applied with the `solve_Neumann` method.
# Note: The number of iterative steps in x-direction is equal to the number of steps in the
# y-direction and taking the argument `sargs`. The time steps takes the argument `targs` with
# step size denoted by `dt`.
# Note: dx = dy. Denoted by the parameter `dx`.
# Note: `diffusion` is the diffusion coefficent and is constant.
# Note: `init` is the initial condition
class diffusion_2dims:
    def __init__(self, targs: int, sargs: int, dt: float, dx: float, diffusion: float, init: np.ndarray) -> None:
        self.targs = targs
        self.sargs = sargs
        self.dt = dt
        self.dx = dx
        self.diffusion = diffusion
        self.init = init

    # Solves the 2d linear diffusion equation with Dirichlet boundary conditions
    def solve_Dirichlet(self, boundary: list) -> np.ndarray:
        try:
            if len(boundary) == 4:  
                temp = []
                temp.append(self.init)
                for iargs in range(0, self.targs):

                    # sets up next arrage of zeros in solution list 
                    next_array = np.zeros((self.sargs, self.sargs))
                    temp.append(next_array)

                    # Finite difference method in 2 dimensions
                    temp[iargs + 1] = temp[iargs] + ((self.diffusion*self.dt)/self.dx)*(
                    np.roll(temp[iargs], 1, axis=1) + np.roll(temp[iargs], -1, axis=1) + 
                    np.roll(temp[iargs], 1, axis=0) + np.roll(temp[iargs], -1, axis=0) - 4*temp[iargs]
                    )

                    # Dirichlet boundary conditions
                    temp[iargs][:, 0] = boundary[0]
                    temp[iargs][:, -1] = boundary[1]
                    temp[iargs][0, :] = boundary[2]
                    temp[iargs][-1, :] = boundary[3]
                
                return np.array(temp)
            else:
                print(f'Boundary list which has lenght {len(boundary)} is the wrong size: length should be 4!')
        except:
            print(f'Boundary list which has lenght {len(boundary)} is the wrong size: length should be 4!')

    # Solves the 2d linear diffusion equation with Neumann boundary conditions
    def solve_Neumann(self, boundary_flux: list):
        try:
            if len(boundary_flux) == 4:  
                temp = []
                temp.append(self.init)
                for iargs in range(0, self.targs):

                    # sets up next arrage of zeros in solution list 
                    next_array = np.zeros((self.sargs, self.sargs))
                    temp.append(next_array)

                    # Finite difference method in 2 dimensions
                    temp[iargs + 1] = temp[iargs] + ((self.diffusion*self.dt)/self.dx)*(
                    np.roll(temp[iargs], 1, axis=1) + np.roll(temp[iargs], -1, axis=1) + 
                    np.roll(temp[iargs], 1, axis=0) + np.roll(temp[iargs], -1, axis=0) - 4*temp[iargs]
                    )

                    # Neumann boundary conditions
                    temp[iargs][:, 0] = temp[iargs][:, 1] + boundary_flux[0]*self.dx
                    temp[iargs][:, -1] = temp[iargs][:, -2] + boundary_flux[1]*self.dx
                    temp[iargs][0, :] = temp[iargs][1, :] + boundary_flux[2]*self.dx
                    temp[iargs][-1, :] = temp[iargs][-2, :] + boundary_flux[3]*self.dx
                
                return np.array(temp)
            else:
                print(f'Boundary list which has lenght {len(boundary_flux)} is the wrong size: length should be 4!')
        except:
            print(f'Boundary list which has lenght {len(boundary_flux)} is the wrong size: length should be 4!')


# Uses the finite difference method in 1 dimension to solve the linear diffusion equation.
# Dirichlet boundary conditions are automaticly applied in the `solve_Dirichlet` method.
# Neumann boundary conditions can be applied with the `solve_Neumann` method.
# Note: The number of iterative steps in x-direction take
# the argument `sargs`. The time steps takes the argument `targs` with
# step size denoted by `dt`.
# Note: dx denoted by the parameter `dx`.
# Note: `diffusion` is the diffusion coefficent and is constant.
# Note: `init` is the initial condition
class diffuision_1dims:
    def __init__(self, targs: int, sargs: int, dt: float, dx: float, diffusion: float, init: np.ndarray) -> None:
        self.targs = targs
        self.sargs = sargs
        self.dt = dt
        self.dx = dx
        self.diffusion = diffusion
        self.init = init

    # Solves the 1d linear diffusion equation with Dirichlet boundary conditions
    def solve_Dirichlet(self, boundary: list) -> np.ndarray:
        try:    
            if len(boundary) == 2:  
                temp = []
                temp.append(self.init)
                for iargs in range(0, self.targs):

                    # sets up next arrage of zeros in solution list 
                    next_array = np.zeros(self.sargs)
                    temp.append(next_array)

                    # Finite difference method in 1 dimension
                    temp[iargs + 1] = temp[iargs] + ((self.diffusion*self.dt)/(self.dx**2))*(
                    np.roll(temp[iargs], 1) + np.roll(temp[iargs], -1) - 2*temp[iargs]
                    )

                    # Dirichlet boundary conditions
                    temp[iargs][ 0] = boundary[0]
                    temp[iargs][-1] = boundary[1]
                
                return np.array(temp)
            else:
                print(f'Boundary list which has lenght {len(boundary)} is the wrong size: length should be 2!')
        except:
            print(f'Boundary list which has lenght {len(boundary)} is the wrong size: length should be 2!')

    def solve_Neumann(self, boundary_flux: list):
        try:    
            if len(boundary_flux) == 2:  
                temp = []
                temp.append(self.init)
                for iargs in range(0, self.targs):

                    # sets up next arrage of zeros in solution list 
                    next_array = np.zeros(self.sargs)
                    temp.append(next_array)

                    # Finite difference method in 1 dimension
                    temp[iargs + 1] = temp[iargs] + ((self.diffusion*self.dt)/(self.dx**2))*(
                    np.roll(temp[iargs], 1) + np.roll(temp[iargs], -1) - 2*temp[iargs]
                    )

                    # Dirichlet boundary conditions
                    temp[iargs][0] = temp[iargs][1] + boundary_flux[0]*self.dx
                    temp[iargs][-1] = temp[iargs][-2] + boundary_flux[1]+self.dx
                
                return np.array(temp)
            else:
                print(f'Boundary list which has lenght {len(boundary_flux)} is the wrong size: length should be 2!')
        except:
            print(f'Boundary list which has lenght {len(boundary_flux)} is the wrong size: length should be 2!')