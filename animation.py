import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

# Save command line arguments
arguments = sys.argv[1:]

# Load data file
fname = 'solution.dat'
with open(fname) as f:
    lines = f.readlines()
lines  = [x.strip() for x in lines]

# Inputs
numbers = [int(i) for i in lines[0].split()]
n_cells = numbers[0]
n_times = numbers[1]

# Initialize arrays
data = np.zeros((n_cells, 5, n_times))
times = np.zeros(n_times)

# Load data into data array
for n in range(0,n_times):
    times[n] = float(lines[1 + n*(n_cells+2)])
    for i in range(0,n_cells):
        string = lines[2 + n*(n_cells+2) + i]
        split = string.split()
        numbers = [float(j) for j in split]
        data[i,:,n] = numbers[0:5]

# Choose index to plot
plotting_index = 4

# Make plot, if desired
if 'image' in arguments:
    fig = plt.figure(figsize=(12,5))
    plot_iter = 190
    plt.plot(data[:,0,plot_iter], data[:,plotting_index,plot_iter], color='k', marker='o', linewidth=3, markersize=8)
    plt.xlabel('x (m.)', fontsize=20)
    plt.ylabel('p (Pa)', fontsize=20)
    print(times[plot_iter]);
    fig.savefig("results.pdf", bbox_inches='tight')
    plt.show()
else:

    # Create figure
    fig, ax = plt.subplots(figsize=(15,7))
    line, = ax.plot(data[:,0,0], data[:,plotting_index,0], color='k', marker='o')
    plt.xlabel('x (m.)')
    plt.ylabel('p (Pa)')

    # Define what happens at every frame in animation
    def update(n, plotting_index, data, line):
        line.set_data(data[:,0,n], data[:,plotting_index,n])
        line.axes.axis([0, 1, 5e5, 2e6])
        ax.set_title("Pressure over Distance at t = " + str(times[n]) + " s.")
        return line,

    # Run animation
    ani = animation.FuncAnimation(fig, update, n_times, fargs=[plotting_index, data, line],
                                  interval=10, blit=False)

    # Set up formatting for the movie files
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=1000000)
    #ani.save('results.mp4', writer=writer)

    plt.show()
