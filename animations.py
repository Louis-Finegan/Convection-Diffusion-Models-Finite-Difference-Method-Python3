import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

# Animates the solution curve in 2 dimensions
def animation_2(solution: np.ndarray, X: np.ndarray, Y: np.ndarray, xlab: str, ylab: str, zlab: str, title: str, zlim: tuple, fps: float, frn: int, filename: str):

    def change_plot(frame_num, solution, plot):
        plot[0].remove()
        plot[0] = ax.plot_surface(X, Y, solution[frame_num, :, :], cmap="afmhot_r")

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot = [ax.plot_surface(X, Y, solution[0, :, :], color='0.75', rstride=1, cstride=1)]
    ax.set_zlim(zlim[0], zlim[1])

    # x label
    if xlab is not None:
        ax.set_xlabel(xlab)
    
    # y label
    if ylab is not None:
        ax.set_ylabel(ylab)
    
    # z label
    if zlab is not None:
        ax.set_zlabel(zlab)

    # title
    if title is not None:
        ax.set_title(title)


    ani = animation.FuncAnimation(fig, change_plot, frn, fargs=(solution, plot), interval=1 / fps)
    plt.close()
    ani.save(filename, writer='pillow', fps=fps)


# Animates the solution curve in 1 dimension
def animation_1(solution: np.ndarray, X: np.ndarray, xlab: str, ylab: str, title: str, color: str, xlim_: tuple, ylim_: tuple, fps: float, frn: int, filename: str):
    fig = plt.figure()
    ax = plt.axes(xlim=(xlim_[0], xlim_[1]), ylim=(ylim_[0], ylim_[1]))
    line, = ax.plot(X, solution[0], color=color)

    # animation function.  This is called sequentially
    def animate(i):
        y = solution[i, :]
        global X
        line.set_ydata(y)
        return line,

    # x label
    if xlab is not None:
        plt.xlabel(xlab)
        
    # y label
    if ylab is not None:
        plt.ylabel(ylab)

    # title
    if title is not None:
        plt.title(title)

    # call the animator.
    anim = animation.FuncAnimation(fig, animate, frames=frn, interval=20, blit=True)
    plt.close()
    anim.save(filename, writer='pillow', fps=fps)


# Animates the solution curve in 2 dimensions but as a color plot
def animation_color(solution: np.ndarray, xlab: str, ylab: str, title: str, xlim_: tuple, ylim_: tuple, fps: float, frn: int, filename: str):
    # Create the figure and axis objects
    fig, ax = plt.subplots()

    # Set the axis limits
    ax.set_xlim(xlim_[0], xlim_[1])
    ax.set_ylim(ylim_[0], ylim_[1])

    # Define the plotting function
    def update_plot(i):
        im.set_array(solution[i, :, :])

    # Create the initial plot
    im = ax.imshow(solution[0, :, :], extent=[xlim_[0], xlim_[1], ylim_[0], ylim_[1]], cmap='viridis')

    # x label
    if xlab is not None:
        plt.xlabel(xlab)
        
    # y label
    if ylab is not None:
        plt.ylabel(ylab)

    # title
    if title is not None:
        plt.title(title)

    # Create the animation object
    anim = animation.FuncAnimation(fig, update_plot, frames=frn, interval=1/fps)

    plt.close()
    anim.save(filename, writer='pillow', fps=fps)

# Animates a histogram of `data` with solution curve `solution`.
def animate_histogram(data: np.ndarray, solution: np.ndarray, X: np.ndarray, bins: int, interval: int, xlim: tuple, xlab: str, title: str, color: str, color_curve: str, fps: float, frn: float, filename: str):
    fig, ax = plt.subplots()

    n_cols = data.shape[1]  # Get the number of columns

    if n_cols == 0:
        raise ValueError(f'Number of columns in `data` should be greater than 0. Number of columns given {n_cols}')
    
    if xlim is not None:
        range1, range2 = xlim
    else:
        range1, range2 = data.min(), data.max()

    # Initialize the histogram with the first column of data
    hist, bin_edges = np.histogram(data[:, 0], bins=bins, range=(range1, range2), density=True)
    bar = ax.bar(bin_edges[:-1], hist, width=bin_edges[1] - bin_edges[0], color=color, edgecolor='black', linewidth=1)

    y = solution[0]
    curve, = ax.plot(X, y, color=color_curve)

    def update(frame):
        # Check if the frame is within the valid range of columns
        if frame >= n_cols:
            return

        # Compute the histogram for the current frame's data column
        hist, _ = np.histogram(data[:, frame], bins=bin_edges, density=True)

        # Update the heights of the existing bars in the histogram
        for h, rect in zip(hist, bar):
            rect.set_height(h)

        # Updates each frame of the solution curve
        y = solution[frame, :]
        curve.set_ydata(y)


    # x label
    if xlab is not None:
        plt.xlabel(xlab)

    # title
    if title is not None:
        plt.title(title)
    
    ani = animation.FuncAnimation(fig, update, frames=range(n_cols), interval=interval, repeat=False)
    

    plt.close()
    ani.save(filename, writer='pillow', fps=fps)

