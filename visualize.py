from matplotlib import pyplot as plt
from itertools import zip_longest

step_num = 0
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
print("Data loaded.")
print("Converting to eventplot format...")
steps_actual = [[] for _ in range(77169)]

for s_index, s_value  in enumerate(steps):
    for neuron in s_value:
        steps_actual[neuron].append(s_index)

del steps
# steps = list(map(list, zip_longest(*steps, fillvalue=None)))


print("Plotting...")
plt.eventplot(steps_actual[:2000], colors='black', lineoffsets=1, linelengths=0.3, linewidths=0.3)
plt.show()

