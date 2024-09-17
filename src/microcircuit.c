#include "microcircuit.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
    Common functions
*/

int synapse_number = 0;

// TODO: add a license here
double randn(double mu, double sigma)
{
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;

    if (call == 1)
    {
        call = !call;
        return (mu + sigma * (double)X2);
    }

    do
    {
        U1 = -1 + ((double)rand() / RAND_MAX) * 2;
        U2 = -1 + ((double)rand() / RAND_MAX) * 2;
        W = pow(U1, 2) + pow(U2, 2);
    } while (W >= 1 || W == 0);

    mult = sqrt((-2 * log(W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mu + sigma * (double)X1);
}

float generate_delay(MicrocircuitLayer layer)
{
    return (float)randn(delta_i[layer][0], delta_i[layer][1]);
}

float generate_synaptic_amp(MicrocircuitLayer layer)
{
    return (float)randn(w_i[layer][0], w_i[layer][1]);
}

float generate_initial_potential(MicrocircuitLayer layer)
{
    return (float)randn(u_init[layer][0], u_init[layer][1]);
}

float generate_uniform_probability()
{
    return (float)rand() / RAND_MAX;
}

float generate_poisson_probability(float lambda, float k)
{
    return (pow(lambda, k) * exp(-lambda)) / tgamma(k + 1);
}

int generate_thalamic_spikes(MicrocircuitLayer layer)
{
    int spikes = 0;
    for (int i = 0; i < thalamic_sizes[layer]; ++i)
    {
        if (generate_uniform_probability() < F_TH * TIMESTEP)
            spikes++;
    }
    return spikes;
}

/*
    Neuron
*/

void create_neuron(LIFNeuron *neuron, MicrocircuitLayer layer)
{
    neuron->membrane = generate_initial_potential(layer);
    neuron->synaptic_amp = generate_synaptic_amp(layer);
    neuron->delay = generate_delay(layer);
    neuron->spike = 0;
    neuron->refractory = 0;
    neuron->synapses = (LIFSynapse *)malloc(INITIAL_SYNAPSE_NUMBER * sizeof(LIFSynapse)); // this is only for the beginning (20000 synapses will NEVER be made)
}

void generate_presynaptic_current(LIFNetwork *network, LIFNeuron *neuron, MicrocircuitLayer layer)
{
    LIFNeuron *neuron_holder;
    neuron->presynaptic_current *= P_11;
    // presynaptic neurons
    for (int i = 0; i < neuron->synapse_count; ++i)
    {
        neuron_holder = network->layers[neuron->synapses[i].pre_layer] + neuron->synapses[i].pre_index;
        if (neuron_holder->spike)
        {
            neuron->presynaptic_current += neuron_holder->synaptic_amp;
        }
    }

    // thalamic inputs
    neuron->presynaptic_current += F_TH * thalamic_sizes[layer] * W_EXT * TAU_SYN; // for now not this: generate_thalamic_spikes(layer) * W_F * W_EXT;
}

void update_neuron(LIFNetwork *network, LIFNeuron *neuron, MicrocircuitLayer layer)
{
    // assume spike = 0
    neuron->spike = 0;

    // update current
    generate_presynaptic_current(network, neuron, layer);

    // update refractory
    if (neuron->membrane >= U_THR)
    {
        neuron->spike = 1;
        neuron->refractory = TAU_REF;
    }
    else if (neuron->refractory > 0)
        neuron->refractory--;
    else
        neuron->refractory = 0;

    // update membrane
    if (neuron->refractory > 0)
        neuron->membrane = U_REST;
    else
        neuron->membrane = P_22 * neuron->membrane + P_21 * neuron->presynaptic_current;
}

/*
    Synapse
*/

void generate_synapse(LIFNeuron *neuron, MicrocircuitLayer pre_layer, int pre_index)
{
    neuron->synapses[neuron->synapse_count].pre_layer = pre_layer;
    neuron->synapses[neuron->synapse_count].pre_index = pre_index;
    neuron->synapse_count++;
}

void create_synapses(LIFNetwork *network)
{
    // int counter = 0;
    // TODO: those loops can be optimized (indexing)
    for (int post_layer = 0; post_layer < LAYER_NUMBER; ++post_layer)
    {
        for (int post_neuron = 0; post_neuron < pop_sizes[post_layer]; ++post_neuron)
        {
            printf("\rGenerating synapses for neuron %d in layer %d", post_neuron, post_layer);
            for (int pre_layer = 0; pre_layer < LAYER_NUMBER; ++pre_layer)
            {
                for (int pre_neuron = 0; pre_neuron < pop_sizes[pre_layer]; ++pre_neuron)
                {
                    if (generate_uniform_probability() <= synaptic_probability[post_layer][pre_layer])
                    {
                        generate_synapse(network->layers[post_layer] + post_neuron, pre_layer, pre_neuron);
                        synapse_number++;
                    }
                }
            }
            // reallocating after every neuron has their connections made to save space
            (network->layers[post_layer] + post_neuron)->synapses = (LIFSynapse *)realloc((network->layers[post_layer] + post_neuron)->synapses, (network->layers[post_layer] + post_neuron)->synapse_count * sizeof(LIFSynapse));
        }
        printf("\n");
    }
}

/*
    Network
*/

void initialize_network(LIFNetwork *network)
{
    network->layers = (LIFNeuron **)malloc(LAYER_NUMBER * sizeof(LIFNeuron *));
    for (int i = 0; i < LAYER_NUMBER; ++i)
    {
        printf("Initializing layer %d: ", i);
        *(network->layers + i) = (LIFNeuron *)malloc(pop_sizes[i] * sizeof(LIFNeuron));
        for (int j = 0; j < pop_sizes[i]; ++j)
        {
            create_neuron((network->layers[i]) + j, i);
        }
        printf("DONE\n");
    }
    create_synapses(network);
    printf("Synapses generated: %d\n", synapse_number);
}

void deinitialize_network(LIFNetwork *network)
{
    for (int i = 0; i < LAYER_NUMBER; ++i)
    {
        for (int j = 0; j < (sizeof(network->layers) / sizeof(LIFNeuron)); ++j)
        {
            free((network->layers[i] + j)->synapses);
        }
        free(network->layers[i]);
    }
    free(network->layers);
}

void update_network(LIFNetwork *network)
{
    for (int i = 0; i < LAYER_NUMBER; ++i)
    {
        for (int j = 0; j < pop_sizes[i]; ++j)
        {
            update_neuron(network, network->layers[i] + j, i);
        }
    }
}

void save_spikes(LIFNetwork *network, int step)
{
    FILE *f;
    f = fopen("spikes.txt", "a+");
    fprintf(f, "step:%d\n", step);
    for (int i = 0; i < LAYER_NUMBER; ++i)
    {
        fprintf(f, "layer:%d\n", i);
        for (int j = 0; j < pop_sizes[i]; ++j)
        {
            fprintf(f, "%d,", (network->layers[i] + j)->spike);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}