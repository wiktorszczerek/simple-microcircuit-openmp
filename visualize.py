from tabnanny import verbose
from matplotlib import pyplot as plt
step_num = 0
TIMESTEP = 0.0001
step = []
steps = []

print("Loading spikes from spikes.txt...")
with open("./spikes.txt", "r") as file:
    for line in file:
        if 'step' in line:
            if step:
                steps.append(step)
            step_num = step_num + 1
        else:
            step = [s.strip('\n') for s in line.split(',')]
            step = [s for s in step if s]
            step = [int(i) for i in step]

del step
# print("Data loaded.")
# print("Plotting the spike raster")
# colors1 = [f'C{i}' for i in range(2)]
# plt.eventplot(steps, colors='black', lineoffsets=1, linelengths=0.2, linewidths=0.5, orientation='vertical', antialiased=True, snap=False)
# plt.axhline(y=20683, color='r', linestyle='--')
# plt.axhline(y=26517, color='g', linestyle='--')
# plt.axhline(y=48432, color='b', linestyle='--')
# plt.axhline(y=53911, color='c', linestyle='--')
# plt.axhline(y=58761, color='m', linestyle='--')
# plt.axhline(y=59826, color='y', linestyle='--')
# plt.axhline(y=74221, color='k', linestyle='--')
# plt.ylim(bottom=0, top=77169)
# plt.xlim(left=0, right=step_num)
# plt.show()


print("Extracting additional data")
neurons = [[] for _ in range(77169)]

for s_index, s_value  in enumerate(steps):
    for neuron in s_value:
        neurons[neuron].append(s_index)

del(steps)

freqs = [[], [], [], [], [], [], [], []]
print(len(freqs))
for index,neuron in enumerate(neurons):
    if(index < 20683):
        freqs[0].append(len(neuron) / (TIMESTEP * step_num))
    elif(index < 26517):
        freqs[1].append(len(neuron) / (TIMESTEP * step_num))
    elif(index < 48432):
        freqs[2].append(len(neuron) / (TIMESTEP * step_num))
    elif(index < 53911):
        freqs[3].append(len(neuron) / (TIMESTEP * step_num))
    elif(index < 58761):
        freqs[4].append(len(neuron) / (TIMESTEP * step_num))
    elif(index < 59826):
        freqs[5].append(len(neuron) / (TIMESTEP * step_num))
    elif(index < 74221):
        freqs[6].append(len(neuron) / (TIMESTEP * step_num))
    else:
        freqs[7].append(len(neuron) / (TIMESTEP * step_num))

print("Plotting additional data")
binwidth = 10
plt.plot(freqs[1], '.')
plt.show()
# del steps
# steps = list(map(list, zip_longest(*steps, fillvalue=None)))



