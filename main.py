import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib.animation as animation
import diffusion as d

def animation_(solution, X, Y, fps, frn) -> None:
    plt.rcParams["figure.figsize"] = [10, 10]
    plt.rcParams["figure.autolayout"] = True

    def change_plot(frame_num, solution, plot):
        plot[0].remove()
        plot[0] = ax.plot_surface(X, Y, solution[frame_num, :, :], cmap="afmhot_r")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot = [ax.plot_surface(X, Y, solution[0, :, :], color='0.75', rstride=100, cstride=100)]
    ax.set_zlim(-5, 10)
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_zlabel('Temperature (10*K)')
    ani = animation.FuncAnimation(fig, change_plot, frn, fargs=(solution, plot), interval=1 / fps)
    plt.show()
    ani.save('animation_convection_diffusion.gif', writer='pillow', fps=fps)

if __name__ == '__main__':
    print('-------------------- \nSCRIPT HAS STARTED \n--------------------')
    dx = 0.01
    dt = 0.01
    x = np.arange(-5, 5, dx)
    y = np.arange(-5, 5 ,dx)
    X, Y = np.meshgrid(x, y)
    init = 10*(1/np.sqrt(0.01*2*np.pi))*np.exp(-(1/2)*((X**2 + Y**2)/0.01))
    heat_array = d.diffusion_2dims(100, 100, dt, dx, 0.1, init)
    print('CLASS INITIALIZED SUCCESSFULLY \n--------------------')
    solution = heat_array.solve()
    print('SOLUTION FOUND')
    animation_(solution, X, Y, 100, 100)
    print('PLOT RAN SUCCESSFULLY \n--------------------')
    print('SCRIPT FINISHED \nSCRIPT RAN SUCCESSFULLY \n--------------------')


