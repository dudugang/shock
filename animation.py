import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Load data file
fname = 'solution.dat'
with open(fname) as f:
    lines = f.readlines()
lines  = [x.strip() for x in lines]

# Inputs
numbers = [int(i) for i in lines[0].split()]
n_cells = numbers[0];
n_times = numbers[1];

# Initialize arrays
data = np.zeros((n_cells, 8, n_times))
times = np.zeros(n_times)

# Load data into data array
for n in range(0,n_times):
    times[n] = float(lines[1 + n*(n_cells+2)]);
    for i in range(0,n_cells):
        string = lines[2 + n*(n_cells+2) + i]
        split = string.split()
        numbers = [float(j) for j in split]
        data[i,:,n] = numbers[0:8]

# Create figure
fig, ax = plt.subplots()
line, = ax.plot(data[:,0,0], data[:,7,0], color='k', marker='o')
plt.xlabel('x (m.)')
plt.ylabel('p (Pa)')

# Define what happens at every frame in animation
def update(n, data, line):
    line.set_data(data[:,0,n], data[:,7,n])
    line.axes.axis([0, 10, 0, 200])
    return line,

# Run animation
ani = animation.FuncAnimation(fig, update, n_times, fargs=[data, line],
                              interval=25, blit=True)

# Save gif
# ani.save('test.gif')
plt.show()
