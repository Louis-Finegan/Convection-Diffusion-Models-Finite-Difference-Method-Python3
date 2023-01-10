

import numpy as np


class diffusion_2dims:
    def __init__(self, targs, sargs, dt, dx, diffusion, init) -> None:
        self.targs = targs
        self.sargs = sargs
        self.dt = dt
        self.dx = dx
        self.diffusion = diffusion
        self.init = init

    def solve(self) -> np.ndarray:
        temp = []
        temp.append(self.init)
        for iargs in range(0, self.targs):
            next_array = np.zeros((self.sargs, self.sargs))
            temp.append(next_array)
            temp[iargs + 1] = temp[iargs] + ((self.diffusion*self.dt)/self.dx)*(
            np.roll(temp[iargs], 1, axis=1) + np.roll(temp[iargs], -1, axis=1) + 
            np.roll(temp[iargs], 1, axis=0) + np.roll(temp[iargs], -1, axis=0) - 4*temp[iargs]
            )

            temp[iargs][:, 0] = 0
            temp[iargs][:, -1] = 0
            temp[iargs][0, :] = 0 
            temp[iargs][-1, :] = 0
        return np.array(temp)

    