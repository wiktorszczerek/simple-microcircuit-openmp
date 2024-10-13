#include "microcircuit.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_randist.h>

/**************************
    Common functions
***************************/
/*
    GSL PRNG.
*/
gsl_rng *gsl_gen;

double generate_delay(MicrocircuitLayer layer)
{
    double res;
    while ((res = gsl_ran_gaussian(gsl_gen, delta_i[layer][1]) + delta_i[layer][0]) <= 0.0)
    {
    };
    return res;
}

double generate_synaptic_amp(MicrocircuitLayer layer)
{
    double res;
    while ((res = gsl_ran_gaussian(gsl_gen, w_i[layer][1]) + w_i[layer][0]) == 0.0)
    {
    };
    return res;
}

double generate_initial_potential(MicrocircuitLayer layer)
{
    double res;
    while ((res = gsl_ran_gaussian(gsl_gen, u_init[layer][1]) + u_init[layer][0]) >= 0.0)
    {
    };
    return res;
}

uint32_t get_pseudorandom_int(uint32_t start, uint32_t stop)
{
    return (uint32_t)gsl_ran_flat(gsl_gen, (double)start, (double)stop);
}

uint32_t generate_thalamic_spikes(MicrocircuitLayer layer)
{
    return gsl_ran_poisson(gsl_gen, TIMESTEP * F_TH * thalamic_sizes[layer]);
}

void random_test()
{
    gsl_gen = gsl_rng_alloc(gsl_rng_mt19937); // allocating and init MT PRNG (seed 0)
    gsl_rng_set(gsl_gen, 0);
    uint32_t test_size = 100000;
    double test[test_size];
    double avg = 0.0;
    double variance = 0.0;
    double std = 0.0;

    for (uint8_t j = 0; j < 8; ++j)
    {
        avg = 0.0;
        for (uint32_t i = 0; i < test_size; ++i)
        {
            avg += (double)generate_thalamic_spikes(j);
        }
        printf("POISSON LAY%u: %lf\n", j, avg / test_size);
        fflush(stdout);
    }
    avg = 0.0;

    for (uint32_t i = 0; i < test_size; ++i)
    {
        test[i] = generate_synaptic_amp(1);
        avg += test[i];
    }
    avg /= test_size;
    for (uint32_t i = 0; i < test_size; ++i)
    {
        variance += pow(test[i] - avg, 2);
    }
    variance /= test_size;
    std = sqrt(variance);
    printf("Randn test: mu = %lf; sigma2 = %lf; sigma = %lf\n", avg, variance, std);
    fflush(stdout);
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

/**
 * Saving data
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

void save_spiking_rates(LIFNetwork *network)
{
    FILE *f;
    f = fopen("additional_data.txt", "w+");
    uint32_t spike_nums[8] = {0};
    for (uint32_t i = 0; i < NEURON_NUMBER; ++i)
    {
        spike_nums[network->neurons[i].layer] += network->neurons[i].spike_number;
        fprintf(f, "Neuron %u: membrane: %f, delay: %u, syn: %f\n", i, network->neurons[i].membrane, network->neurons[i].delay, network->neurons[i].synaptic_amp);
    }
    for (uint8_t i = 0; i < 8; ++i)
    {
        fprintf(f, "Layer %d: %u\n", i, spike_nums[i]);
    }
    fclose(f);
}

/**************************
        Neuron
***************************/

void create_neuron(LIFNetwork *network, MicrocircuitLayer layer, uint32_t neuron_index)
{
    double delay;
    LIFNeuron *neuron = &(network->neurons[pop_starts[layer] + neuron_index]);
    neuron->layer = layer;
    neuron->membrane = generate_initial_potential(layer);
    neuron->synaptic_amp = generate_synaptic_amp(layer);
    delay = ceil(10 * generate_delay(layer));
    neuron->delay = delay;
    neuron->spike = 0;
    neuron->refractory = 0;
    neuron->synapse_count = 0;
    neuron->spike_timestamps[0] = 0;
    neuron->spike_timestamps[1] = 0;
    neuron->total_current = 0.0;
    neuron->presynaptic_current = 0.0;
    neuron->spike_timestamp_flag = 0;
    neuron->spike_number = 0;
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
        Update
***************************/
void update_neuron(LIFNeuron *neuron, uint32_t current_timestep, uint32_t thalamic_spikes)
{

    neuron->total_current *= P_11;
    neuron->total_current += neuron->presynaptic_current * W_F; // sum = spike_ij * weight_ij [V] [F/s] * [V] = [A]
                                                                // update the total current with thalamics
                                                                // decay rate [1] - this only applies to the "old" spikes. Thalamics are from THIS step
    neuron->total_current += thalamic_spikes * W_EXT * W_F;     // thalamic_spikes * W_EXT * W_F; // [F/s] * [V] = [A]

    // update membrane
    if (neuron->refractory > 0)
        neuron->membrane = U_REST;
    else
        neuron->membrane = P_22 * neuron->membrane + P_21 * neuron->total_current;

    // assume spike = 0 and update refractory
    neuron->spike = 0;
    if (neuron->membrane >= U_THR)
    {
        if (!current_timestep)
            neuron->spike_timestamps[0] = current_timestep + neuron->delay;
        else
            neuron->spike_timestamps[1] = current_timestep + neuron->delay;
        neuron->spike = 1;
        neuron->spike_number++;
        neuron->refractory = TAU_REF;
    }
    else if (neuron->refractory > 0)
        neuron->refractory--;
}

void update_network(LIFNetwork *network, uint32_t current_timestep)
{
    uint32_t neuron, pre_neuron;
    LIFNeuron *current_neuron, *pre_neuron_ptr;
    network->current_timestep = current_timestep;

#ifdef MULTIPROCESSING
    uint32_t thalamic_spikes[8];
    uint32_t thread_id;

// thalamics here are from t-1
#pragma omp parallel shared(thalamic_spikes, current_timestep) num_threads(8) private(thread_id)
    {
        thread_id = omp_get_thread_num();
        gsl_rng_set(gsl_gen, thread_id);
        thalamic_spikes[thread_id] = generate_thalamic_spikes(thread_id);
#pragma omp barrier
    }

#ifdef MULTIPROCESSING
// presynptics are from t - 1 as well
#pragma omp parallel for private(pre_neuron, neuron, pre_neuron_ptr, current_neuron) shared(network)
#endif
    for (neuron = 0; neuron < NEURON_NUMBER; ++neuron)
    {
        current_neuron = &(network->neurons[neuron]);
        current_neuron->presynaptic_current = 0.0; // this needs to be zeroed!
        for (pre_neuron = 0; pre_neuron < network->neurons[neuron].synapse_count; ++pre_neuron)
        {
            pre_neuron_ptr = current_neuron->presynaptic_neurons[pre_neuron];
            if ((pre_neuron_ptr->spike_timestamps[0] - network->current_timestep) == 0 && pre_neuron_ptr->spike_timestamps[0]) // exactly this
            {
                current_neuron->presynaptic_current += pre_neuron_ptr->synaptic_amp;
                pre_neuron_ptr->spike_timestamp_flag = 1; // doesn't need to be explicitly atomic
            }
        }
    }

#pragma omp parallel for private(neuron) shared(network, thalamic_spikes)
#endif
    for (neuron = 0; neuron < NEURON_NUMBER; ++neuron)
    {
        if (network->neurons[neuron].spike_timestamp_flag)
        {
            network->neurons[neuron].spike_timestamps[0] = network->neurons[neuron].spike_timestamps[1];
            network->neurons[neuron].spike_timestamp_flag = 0;
        }
        update_neuron(&(network->neurons[neuron]), network->current_timestep, thalamic_spikes[network->neurons[neuron].layer]);
    }
}

/*
    Init and deinit
*/
void initialize_network(LIFNetwork *network)
{
    printf("Reserving space for network...");
    network->neurons = (LIFNeuron *)malloc(NEURON_NUMBER * sizeof(LIFNeuron));
    printf("DONE\n");
    printf("Creating neurons...");
    fflush(stdout);
    uint8_t pre_layer;
    uint32_t pre_neuron;

    gsl_gen = gsl_rng_alloc(gsl_rng_mt19937); // allocating and init MT PRNG (seed 0)
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