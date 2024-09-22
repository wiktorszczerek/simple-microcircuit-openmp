import json
from os import path
from sys import argv
import numpy as np
from tqdm import tqdm
import multiprocessing

"""
Structure for config JSON:
{
    "scale": 1,
    "approximate_thalamic": true,
    "pop_sizes": [
        20683, 5834...
    ],
    ...
}
"""

class MicrocircuitGenerator:
    config_filename = 'generator_config.json'
    output_filename = 'network.json'
    data = {}
    network = []

    def __init__(self, config_file_path="./src"):
        with open(path.join(config_file_path, self.config_filename)) as file:
            self.data = json.load(file)

    def generate_layer(self, population_data) -> list:
        """
        neurons is a list of lists (neurons) -> [neurons] = [ [u_init, w_i, delta_i], ... ]
        K will be omitted because I'll model it as a constant value following D&P
        """
        neurons = []
        for neuron in range(0, population_data["size"]):
            neurons.append([
                float(np.random.normal(population_data["u_init"][0], population_data["u_init"][1], 1)[0]),
                float(np.random.normal(population_data["w_i"][0], population_data["w_i"][1], 1)[0]),
                float(np.random.normal(population_data["delta_i"][0], population_data["delta_i"][1], 1)[0]),
            ])
        return neurons
    
    def generate_network(self) -> None:
        """
        Generate the microcircuit in two parts: neuron list and synpase list.
        Should be easier to save it in files and read in C this way.
        """
        if not self.data:
            raise Exception("Data not loaded properly.")
        
        for population_data_key, population_data in self.data["populations"].items():
            self.network.append(self.generate_layer(population_data))

        conns = []
        for post_layer_index, (post_layer_key, post_layer) in enumerate(self.data["populations"].items()):
            for post_neuron_index in tqdm(range(0, post_layer["size"])):
                for pre_layer_index, (pre_layer_key, pre_layer) in enumerate(self.data["populations"].items()):
                    for pre_neuron_index in range(0, pre_layer["size"]):
                        if(np.random.random_sample(1) < self.data["probability_array"][pre_layer_index][post_layer_index]):
                            conns.append([pre_layer_index, pre_neuron_index, post_layer_index, post_neuron_index])

    # def generate_connections(self, pre_layer_index) -> list:
    #     conns = []
    #     for post_pop_index, post_pop in enumerate(self["populations"]):
    #         for post_neuron_index in range(0, post_pop["size"]):
    #             if(np.random.random_sample(1) < self.data["probability_array"][pre_layer_index][post_pop_index]):
    #                 conns.append(post_neuron_index)
    #     return conns

    # def generate_network(self):
    #     if not self.data:
    #         raise Exception("Data not loaded properly.")

    #     population_dest_index = 0
    #     for population_dest_key, population_dest in self.data["populations"].items():
    #         layer = []
    #         for neuron_dest in range(0, population_dest["size"]):
    #             synapses = {}
    #             neuron = {
    #                 "u_init": float(np.random.normal(population_dest["u_init"][0], population_dest["u_init"][1], 1)[0]),
    #                 "w_i": float(np.random.normal(population_dest["w_i"][0], population_dest["w_i"][1], 1)[0]),
    #                 "delta_i": float(np.random.normal(population_dest["delta_i"][0], population_dest["delta_i"][1], 1)[0]),
    #                 "synapses": {}
    #             }
                # population_src_index = 0
                # for population_src_key, population_src in self.data["populations"].items():
                #     layer_synapses = []
                #     for neuron_src in range(0, population_src["size"]):
                #         if(np.random.random_sample(1) < self.data["probability_array"][population_src_index][population_dest_index]):
                #             layer_synapses.append(neuron_src)
                #     synapses[population_src_key] = layer_synapses
                #     population_src_index += 1
            #     neuron["synapses"] = synapses
            #     layer.append(neuron)
            # population_dest_index += 1
            # self.network[population_dest_key] = layer
    
    # def save_network(self):
    #     with open(path.join(config_file_path, self.config_filename)) as file:
    #         self.data = json.load(file)

def main():
    gen = MicrocircuitGenerator()
    gen.generate_network()


if __name__ == "__main__":
    main()
