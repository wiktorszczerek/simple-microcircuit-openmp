#include "microcircuit.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <assert.h>

/**************************
    Common functions
***************************/

/*
    Implementation of randn() published on: https://kcru.lawsonresearch.ca/research/srk/normalDBN_random.html.
    All rights go to their respective owners.
*/
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

uint32_t get_pseudorandom_int(uint32_t start, uint32_t stop)
{
    return rand() % (stop + 1 - start) + start;
}

double generate_uniform_probability()
{
    return ((double)rand() / (double)RAND_MAX);
}

float generate_poisson_probability(float lambda, float k)
{
    return (pow(lambda, k) * exp(-lambda)) / tgamma(k + 1);
}

uint32_t generate_thalamic_spikes(MicrocircuitLayer layer)
{
    uint32_t spikes = 0;
    for (uint32_t i = 0; i < thalamic_sizes[layer]; ++i)
    {
        if (generate_uniform_probability() < F_TH * TIMESTEP)
            spikes++;
    }
    return spikes;
}

uint32_t max_synaptic_number_per_layer(MicrocircuitLayer layer)
{
    uint32_t num = 0;
    for (uint32_t i = 0; i < LAYER_NUMBER; ++i)
    {
        num += max_synapses_per_layer[layer][i];
    }
    return num;
}

uint32_t synaptic_number_check()
{
    printf("==============================================================================================================\n");
    printf("Probability array check. Checking for postsynaptic connections per neuron\n");
    uint32_t synapse_number = 0, maximum = 0, synapses_per_population = 0;
    for (uint32_t i = 0; i < LAYER_NUMBER; ++i)
    {
        for (uint32_t j = 0; j < LAYER_NUMBER; ++j)
        {
            synapses_per_population = (uint32_t)ceil((double)synaptic_probability[i][j] * pop_sizes[j] * pop_sizes[i]);
            if (synapses_per_population > maximum)
                maximum = synapses_per_population;
            synapse_number += synapses_per_population;
            printf("%12d ", synapses_per_population);
        }
        printf("\n");
    }
    printf("==============================================================================================================\n");
    return maximum;
}

uint32_t get_linux_random()
{
    FILE *file = fopen("/dev/urandom", "r");
    uint32_t i = 0;
    if (!fread(&i, sizeof(i), 1, file))
        printf("ERROR: Couldn't read the Linux random. Defaulting to 0.\n");
    fclose(file);
    return i;
}

uint32_t get_neuron_index(MicrocircuitLayer layer, uint32_t neuron_index)
{
    uint32_t ret = 0;
    for (uint32_t i = 0; i < layer; ++i)
    {
        ret += pop_sizes[i];
    }
    ret += neuron_index;
    return ret;
}

/**
 * Saving the absolute indexes [0...77169] of the spiked neuron in the current timestep
 */
void save_spikes(LIFNetwork *network, uint32_t step)
{
    FILE *f;
    f = fopen("spikes.txt", "a+");
    fprintf(f, "step:%d\n", step);
    for (uint32_t i = 0; i < NEURON_NUMBER; ++i)
    {
        if (network->neurons[i].spike)
            fprintf(f, "%d,", i);
    }
    fprintf(f, "\n");
    fclose(f);
}

/**************************
        Neuron
***************************/

void create_neuron(LIFNetwork *network, MicrocircuitLayer layer, uint32_t neuron_index)
{
    LIFNeuron *neuron = &(network->neurons[pop_starts[layer] + neuron_index]);
    neuron->layer = layer;
    neuron->membrane = generate_initial_potential(layer);
    neuron->synaptic_amp = generate_synaptic_amp(layer);
    neuron->delay = generate_delay(layer);
    neuron->spike = 0;
    neuron->refractory = 0;
    neuron->synapse_count = 0;
}

void update_neuron(LIFNeuron *neuron)
{
    // assume spike = 0
    neuron->spike = 0;

    // thalamic inputs approximation
    neuron->presynaptic_current += F_TH * thalamic_sizes[neuron->layer] * W_EXT * TAU_SYN;

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

    neuron->presynaptic_current *= P_11;
}

/**************************
        Synapse
***************************/

void create_synapse_pairs(LIFNetwork *network)
{
    network->synapse_count = 0;
    for (uint8_t i = 0; i < LAYER_NUMBER; ++i)
    {
        for (uint8_t j = 0; j < LAYER_NUMBER; ++j)
        {
            network->synapse_count += max_synapses_per_layer[i][j];
        }
    }
    LIFConnection *synapses = malloc(network->synapse_count * sizeof(LIFConnection));

    printf("\n");
    // this could be parallelized but there's no point - it only runs once
    uint32_t synapse_index = 0;
    for (uint8_t pre_layer = 0; pre_layer < LAYER_NUMBER; ++pre_layer)
    {
        for (uint8_t post_layer = 0; post_layer < LAYER_NUMBER; ++post_layer)
        {
            printf("Generating synapses for pop_%d -> pop_%d...", pre_layer, post_layer);
            if (!max_synapses_per_layer[pre_layer][post_layer])
            {
                printf("FAILED, no synapses between those populations. Skipping.\n");
                continue;
            }
            for (uint32_t sample = 0; sample < max_synapses_per_layer[pre_layer][post_layer]; ++sample)
            {
                synapses[synapse_index].pre_index = get_pseudorandom_int(0, pop_sizes[pre_layer] - 1);
                synapses[synapse_index].post_index = get_pseudorandom_int(0, pop_sizes[post_layer] - 1);
                synapses[synapse_index].pre_layer = pre_layer;
                synapses[synapse_index].post_layer = post_layer;
                network->neurons[pop_starts[post_layer] + synapses[synapse_index].post_index].synapse_count++;
                synapse_index++;
            }
            printf("DONE\n");
        }
        printf("\n");
    }

    LIFNeuron *pre_neuron, *post_neuron;

    for (uint32_t neuron = 0; neuron < NEURON_NUMBER; ++neuron)
    {
        post_neuron = network->neurons + neuron;
        post_neuron->presynaptic_neurons = malloc(post_neuron->synapse_count * sizeof(LIFNeuron *));
    }

    uint32_t *neuron_indexes = calloc(NEURON_NUMBER, sizeof(uint32_t));
    uint32_t pre_index, post_index;
    for (uint32_t synapse = 0; synapse < network->synapse_count; ++synapse)
    {
        pre_index = pop_starts[synapses[synapse].pre_layer] + synapses[synapse].pre_index;
        post_index = pop_starts[synapses[synapse].post_layer] + synapses[synapse].post_index;
        pre_neuron = &(network->neurons[pre_index]);
        post_neuron = &(network->neurons[post_index]);
        post_neuron->presynaptic_neurons[neuron_indexes[post_index]] = pre_neuron;
        neuron_indexes[post_index]++;
    }
    free(neuron_indexes);
    free(synapses);
}

/**************************
        Network
***************************/

void initialize_network(LIFNetwork *network)
{

#ifdef MULTIPROCESSING
    omp_set_num_threads(NUM_THREADS);
#endif

    printf("Reserving space for network...");
    network->neurons = (LIFNeuron *)malloc(NEURON_NUMBER * sizeof(LIFNeuron));
    printf("DONE\n");
    printf("Creating neurons...");
    fflush(stdout);
    uint8_t pre_layer;
    uint32_t pre_neuron;

#ifdef MULTIPROCESSING
#pragma omp parallel shared(network) num_threads(LAYER_NUMBER)
    srand(0);
#pragma omp for private(pre_layer, pre_neuron)
#endif
    for (pre_layer = 0; pre_layer < LAYER_NUMBER; ++pre_layer)
    {
        for (pre_neuron = 0; pre_neuron < pop_sizes[pre_layer]; ++pre_neuron)
        {
            create_neuron(network, pre_layer, pre_neuron);
        }
    }
    printf("DONE\n");
    printf("Generating neuron pairs");
    create_synapse_pairs(network);
    printf("Generating neuron pairs complete\n");
}

void deinitialize_network(LIFNetwork *network)
{
    printf("Releasing resources...");
    for (uint32_t i = 0; i < NEURON_NUMBER; ++i)
    {
        free(network->neurons[i].presynaptic_neurons);
    }
    free(network->neurons);
    printf("DONE\n");
}

void update_network(LIFNetwork *network)
{
    uint32_t neuron, pre_neuron;
    LIFNeuron *current_neuron, *pre_neuron_ptr;
#ifdef MULTIPROCESSING
#pragma omp parallel for private(pre_neuron, neuron) shared(network)
#endif
    for (neuron = 0; neuron < NEURON_NUMBER; ++neuron)
    {
        current_neuron = &(network->neurons[neuron]);
        for (pre_neuron = 0; pre_neuron < network->neurons[neuron].synapse_count; ++pre_neuron)
        {
            pre_neuron_ptr = current_neuron->presynaptic_neurons[pre_neuron];
            if (pre_neuron_ptr->spike)
                current_neuron->presynaptic_current += pre_neuron_ptr->synaptic_amp;
        }
    }
#ifdef MULTIPROCESSING
#pragma omp barrier
#pragma omp parallel for private(neuron) shared(network)
#endif
    for (neuron = 0; neuron < NEURON_NUMBER; ++neuron)
    {
        update_neuron(&(network->neurons[neuron]));
    }
}
