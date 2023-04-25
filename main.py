import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
import diffuconpy.diffusion as d
import diffuconpy.convection as c
import diffuconpy.convection_diffusion as dc
import time

# Animates the diffusion process in 2 dimensions
def animation_(solution, X, Y, fps, frn, filename):
    plt.rcParams["figure.figsize"] = [10, 10]
    plt.rcParams["figure.autolayout"] = True

    def change_plot(frame_num, solution, plot):
        plot[0].remove()
        plot[0] = ax.plot_surface(X, Y, solution[frame_num, :, :], cmap="afmhot_r")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot = [ax.plot_surface(X, Y, solution[0, :, :], color='0.75', rstride=1, cstride=1)]
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
def animation_color(solution, fps, frn, filename):
    # Create the figure and axis objects
    fig, ax = plt.subplots()

    # Set the axis limits
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)

    # Define the plotting function
    def update_plot(i):
        im.set_array(solution[i, :, :])

    # Create the initial plot
    im = ax.imshow(solution[0, :, :], extent=[-5, 5, -5, 5], cmap='viridis')

    # Create the animation object
    anim = animation.FuncAnimation(fig, update_plot, frames=frn, interval=1/fps)

    plt.close()
    anim.save(filename, writer='pillow', fps=fps)

if __name__ == '__main__':
    # Uses diffusion modules to solve with initial condition `init` and parameters given
    print('-------------------- \nSCRIPT HAS STARTED \n--------------------')
    start_time = time.time()
    dx = 0.1
    dt = 0.1
    x = np.arange(-5, 5, dx)
    y = np.arange(-5, 5 ,dx)
    X, Y = np.meshgrid(x, y)
    init = 4*(1/np.sqrt(0.01*2*np.pi))*np.exp(-(1/2)*((X**2 + Y**2)/0.1))
    #init = np.heaviside(X, 1) * np.heaviside(Y, 1)
    heat_array = d.diffusion_2dims(250, 100, dt, dx, 0.009, init)
    #transport = c.convection_2dims(250, 200, dt, dx, -0.1, init)
    #dc_array = dc.convection_diffusion_2dims(250, 200, dt, dx, 0.009, -0.1, init)
    print('OBJECT INITIALIZED SUCCESSFULLY \n--------------------')
    sol = heat_array.solve_Dirichlet(boundary=[0,0,0,0])
    #sol = transport.solve()
    #sol = dc_array.solve_Dirichlet([0, 0, 0, 0])
    print(time.time() - start_time)
    print(type(sol.solution))
    print(sol)
    if isinstance(sol.solution, np.ndarray):
        print('SOLUTION FOUND')
    else:
        print('--------------------\nSOLUTION NOT FOUND')
    try:   
        animation_color(sol.solution, 60, 250, 'img//animation_diffusion.gif')
        print('PLOT SAVED SUCCESSFULLY \n--------------------')
        print('SCRIPT FINISHED \nSCRIPT RAN SUCCESSFULLY \n--------------------')
    except:
        print('--------------------\nFAILED TO RENDER PLOT\n--------------------\nSCRIPT RAN UNSUCCESSFULLY')

