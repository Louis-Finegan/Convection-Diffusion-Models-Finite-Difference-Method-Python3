import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
import diffuconpy.diffusion as d
import time

# Animates the diffusion process in 2 dimensions
def animation_(solution, X, Y, fps, frn, filename) -> None:
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
    ax.set_zlabel('Temperature (2*K)')
    ani = animation.FuncAnimation(fig, change_plot, frn, fargs=(solution, plot), interval=1 / fps)
    plt.close()
    ani.save(filename, writer='pillow', fps=fps)

# Animates the diffusion process in 1 dimension
def animation_1(solution, X, fps, frn, filename):
        fig = plt.figure()
        ax = plt.axes(xlim=(-5, 5), ylim=(-5, 10)) # left bound -5 and right bound 5
        line, = ax.plot(X, solution[0])

        # animation function.  This is called sequentially
        def animate(i):
            y = solution[i, :]
            global X
            line.set_ydata(y)
            return line,

        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig, animate, frames=frn, interval=20, blit=True)
        plt.close()
        anim.save(filename, writer='pillow', fps=fps)

# Animates the diffusion process in 2 dimensions but as a color plot
def animation_color():
    pass

if __name__ == '__main__':
    # Uses diffusion modules to solve with initial condition `init` and parameters given
    print('-------------------- \nSCRIPT HAS STARTED \n--------------------')
    dx = 0.1
    dt = 0.1
    x = np.arange(-2, 2, dx)
    y = np.arange(-2, 2 ,dx)
    X, Y = np.meshgrid(x, y)
    #init = 4*(1/np.sqrt(0.01*2*np.pi))*np.exp(-(1/2)*((X**2 + Y**2)/0.01))
    init = np.heaviside(X, 1) * np.heaviside(Y, 1)
    heat_array = d.diffusion_2dims(200, 50, dt, dx, 0.1, init)
    print('OBJECT INITIALIZED SUCCESSFULLY \n--------------------')
    solution = heat_array.solve_Dirichlet(boundary=[0,0,0,0])
    print(time.time())
    if isinstance(solution, np.ndarray):
        print('SOLUTION FOUND')
    else:
        print('--------------------\nSOLUTION NOT FOUND')
    try:   
        animation_(solution, X, Y, 20, 200, 'animation_diffusion.gif')
        print('PLOT SAVED SUCCESSFULLY \n--------------------')
        print('SCRIPT FINISHED \nSCRIPT RAN SUCCESSFULLY \n--------------------')
    except:
        print('--------------------\nFAILED TO RENDER PLOT\n--------------------\nSCRIPT RAN UNSUCCESSFULLY')

