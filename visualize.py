from matplotlib import pyplot as plt
import numpy as np

steps = {}

with open("./spikes.txt", "r") as file:
    layers = {}
    for line in file:
        if 'step' in line:
            if layers:
                steps[step_number] = layers
            step_number = line.split(':')[1].strip()
            continue
        elif 'layer' in line:
            layer_number = line.split(':')[1].strip()
            continue
        else:
            spikes = [s.strip('\n') for s in line.split(',')]
            spikes = [s for s in spikes if s]
            spikes = [int(i) for i in spikes]
            layers[layer_number] = spikes
    steps[step_number] = layers

simulation_spikes = []
# step_spikes = []
# spike_times = []

for step in steps.values():
    step_spikes = []
    for layer in step.values():
        step_spikes = step_spikes + layer
    simulation_spikes.append(step_spikes)



simulation_spikes = np.array(simulation_spikes)

simulation_spikes[0][0] = 1
simulation_spikes[1][0] = 1

for i in range(simulation_spikes.shape[0]):
    simulation_spikes[i] = simulation_spikes[i] * (i+1)

simulation_spikes = simulation_spikes.transpose()

sim_spikes = []

print(simulation_spikes.shape)

for i in range(simulation_spikes.shape[0]):
    sim_spikes.append([j for j in simulation_spikes[i] if j != 0])



# for step 

# derp = np.random.gamma(4, size=[60, 50])

f = plt.figure()
plt.eventplot(sim_spikes, colors='black', lineoffsets=1, linelengths=1)
plt.show()
input()
f.clear()
plt.close(f)

