from tabnanny import verbose
from matplotlib import pyplot as plt
import pandas as pd
step_num = 0
TIMESTEP = 0.0001
step = []
steps = []

print("Loading spikes from spikes.txt...")
with open("./spikes.txt", "r") as file:
    for line in file:
        if 'step' in line:
            # if step:
            steps.append(step)
            step_num = step_num + 1
        else:
            step = [s.strip('\n') for s in line.split(',')]
            step = [s for s in step if s]
            step = [int(i) for i in step]

step_num_half = step_num / 2

del step



print("Extracting additional data")
neurons = [[] for _ in range(77169)]

for s_index, s_value  in enumerate(steps):
    for neuron in s_value:
        neurons[neuron].append(s_index)

del(steps)

print("Data loaded.")
print("Plotting the spike raster")
colors1 = [f'C{i}' for i in range(2)]
plt.eventplot(neurons, colors='black', lineoffsets=1, linelengths=0.75, linewidths=0.5, orientation='horizontal', antialiased=False, snap=False)
plt.axhline(y=20683, color='r', linestyle='--')
plt.axhline(y=26517, color='g', linestyle='--')
plt.axhline(y=48432, color='b', linestyle='--')
plt.axhline(y=53911, color='c', linestyle='--')
plt.axhline(y=58761, color='m', linestyle='--')
plt.axhline(y=59826, color='y', linestyle='--')
plt.axhline(y=74221, color='k', linestyle='--')
plt.ylim(bottom=0, top=77169)
plt.xlim(left=0, right=step_num)
plt.show()



freqs = [[], [], [], [], [], [], [], []]
print(len(freqs))
for index,neuron in enumerate(neurons):
    if(index < 20683):
        freqs[0].append(len(neuron)/(TIMESTEP * step_num))
    elif(index < 26517):
        freqs[1].append(len(neuron)/(TIMESTEP * step_num))
    elif(index < 48432):
        freqs[2].append(len(neuron)/(TIMESTEP * step_num))
    elif(index < 53911):
        freqs[3].append(len(neuron)/(TIMESTEP * step_num))
    elif(index < 58761):
        freqs[4].append(len(neuron)/(TIMESTEP * step_num))
    elif(index < 59826):
        freqs[5].append(len(neuron)/(TIMESTEP * step_num))
    elif(index < 74221):
        freqs[6].append(len(neuron)/(TIMESTEP * step_num))
    else:
        freqs[7].append(len(neuron)/(TIMESTEP * step_num))

fig, ax = plt.subplots(4,2)
ax[0,0].plot(freqs[0], '.')
ax[0,1].plot(freqs[1], '.')
ax[1,0].plot(freqs[2], '.')
ax[1,1].plot(freqs[3], '.')
ax[2,0].plot(freqs[4], '.')
ax[2,1].plot(freqs[5], '.')
ax[3,0].plot(freqs[6], '.')
ax[3,1].plot(freqs[7], '.')
plt.show()

# fig, ax = plt.subplots(4,2)
# pd.Series(freqs[0], copy=False).plot.kde(ax=ax[0,0],xlim=0)
# pd.Series(freqs[1], copy=False).plot.kde(ax=ax[0,1],xlim=0)
# pd.Series(freqs[2], copy=False).plot.kde(ax=ax[1,0],xlim=0)
# pd.Series(freqs[3], copy=False).plot.kde(ax=ax[1,1],xlim=0)
# pd.Series(freqs[4], copy=False).plot.kde(ax=ax[2,0],xlim=0)
# pd.Series(freqs[5], copy=False).plot.kde(ax=ax[2,1],xlim=0)
# pd.Series(freqs[6], copy=False).plot.kde(ax=ax[3,0],xlim=0)
# pd.Series(freqs[7], copy=False).plot.kde(ax=ax[3,1],xlim=0)
# plt.xlim(left=0)
# plt.show()



