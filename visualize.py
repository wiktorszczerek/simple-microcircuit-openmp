from matplotlib import pyplot as plt
import numpy as np
from tqdm import tqdm
from pathlib import Path

numpy_array_file = Path('./spikes.npy')


step_num = 0
step = []

steps = {}

if not numpy_array_file.exists():
    print("File spikes.npy not found. Reading spikes from spikes.txt...")

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

    print("Converting spikes...")

    simulation_spikes = []

    for step in steps.values():
        step_spikes = []
        for layer in step.values():
            step_spikes = step_spikes + layer
        simulation_spikes.append(step_spikes)

    del steps, step_spikes, spikes

    simulation_spikes = np.array(simulation_spikes)
    for i in range(simulation_spikes.shape[0]):
        simulation_spikes[i] = simulation_spikes[i] * (i+1)
    simulation_spikes = simulation_spikes.transpose()
    
    print("Data loaded. Saving to spikes.npy...")
    np.save(numpy_array_file, simulation_spikes)
    print("Data saved.")
else:
    print("spikes.npy found. Loading data...")
    simulation_spikes = np.load(numpy_array_file)
    print("Data loaded.")



print("Size of the simulation data: "+str(simulation_spikes.shape))

sim_spikes = []
holder = []
 
print("Converting to eventplot...") # we could do that before but i'd rather save arrays of 0 and 1 rather than empty arrays

for i in tqdm(range(simulation_spikes.shape[0])):
    holder = [j for j in simulation_spikes[i] if j != 0]
    sim_spikes.append(holder)

print("Plotting...")

plt.eventplot(sim_spikes[100:200], colors='black', lineoffsets=1, linelengths=1)
plt.show()

