import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

# Animates the diffusion process in 2 dimensions
def animation_2(solution: np.ndarray, X: np.ndarray, Y: np.ndarray, xlab: str, ylab: str, zlab: str, title: str, zlim: tuple, fps: float, frn: int, filename: str):
    plt.rcParams["figure.figsize"] = [10, 10]
    plt.rcParams["figure.autolayout"] = True

    def change_plot(frame_num, solution, plot):
        plot[0].remove()
        plot[0] = ax.plot_surface(X, Y, solution[frame_num, :, :], cmap="afmhot_r")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot = [ax.plot_surface(X, Y, solution[0, :, :], color='0.75', rstride=1, cstride=1)]
    ax.set_zlim(zlim[0], zlim[1])

    if xlab != None:
        ax.set_xlabel(xlab)
    
    if ylab != None:
        ax.set_ylabel(ylab)
    
    if zlab != None:
        ax.set_zlabel(zlab)

    if title != None:
        ax.set_title(title)


    ani = animation.FuncAnimation(fig, change_plot, frn, fargs=(solution, plot), interval=1 / fps)
    plt.close()
    ani.save(filename, writer='pillow', fps=fps)


# Animates the diffusion process in 1 dimension
def animation_1(solution: np.ndarray, X: np.ndarray, xlab: str, ylab: str, title: str, color: str, xlim_: tuple, ylim_: tuple, fps: float, frn: int, filename: str):
    fig = plt.figure()
    ax = plt.axes(xlim=(xlim_[0], xlim_[1]), ylim=(ylim_[0], ylim_[1]))
    line, = ax.plot(X, solution[0], color = color)

    # animation function.  This is called sequentially
    def animate(i):
        y = solution[i, :]
        global X
        line.set_ydata(y)
        return line,

    if xlab != None:
        plt.xlabel(xlab)
        
    if ylab != None:
        plt.ylabel(ylab)

    if title != None:
        plt.title(title)

    # call the animator.
    anim = animation.FuncAnimation(fig, animate, frames=frn, interval=20, blit=True)
    plt.close()
    anim.save(filename, writer='pillow', fps=fps)


# Animates the diffusion process in 2 dimensions but as a color plot
def animation_color(solution: np.ndarray, fps: float, frn: int, filename: str):
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


def animation_hist(data: np.ndarray, bins: int, color: str, ylim_: float, fps: float, frn: int, filename: str):
    def prepare_animation(bar_container):
        
        def animate(frame_number):
            n, _ = np.histogram(data[:, frame_number], bins)
            for count, rect in zip(n, bar_container.patches):
                rect.set_height(count)
            return bar_container.patches
        return animate

    # Output generated via `matplotlib.animation.Animation.to_jshtml`.
    fig, ax = plt.subplots()
    _, _, bar_container = ax.hist(data, bins, lw=1, ec="white", fc=color, alpha=0.5)
    ax.set_ylim(top=ylim_)  # set safe limit to ensure that all data is visible.

    ani = animation.FuncAnimation(fig, prepare_animation(bar_container), frn, repeat=False, blit=True)

    plt.close()
    ani.save(filename, writer='pillow', fps=fps)