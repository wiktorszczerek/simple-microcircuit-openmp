#include "microcircuit.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_randist.h>
#include <assert.h>

/**************************
    Common functions
***************************/
/*
    GSL PRNG.
*/
#ifdef MULTIPROCESSING
gsl_rng *gsl_gen[16]; // num_threds max
#else
gsl_rng *gsl_gen;
#endif

double generate_delay(MicrocircuitLayer layer)
{
    double res;
#ifdef MULTIPROCESSING
    while ((res = gsl_ran_gaussian(gsl_gen[0], delta_i[layer][1]) + delta_i[layer][0]) <= 0.0)
    {
    };
#else
    while ((res = gsl_ran_gaussian(gsl_gen, delta_i[layer][1]) + delta_i[layer][0]) <= 0.0)
    {
    };
#endif
    return res;
}

double generate_synaptic_amp(MicrocircuitLayer layer)
{
    double res;
#ifdef MULTIPROCESSING
    while ((res = gsl_ran_gaussian(gsl_gen[0], w_i[layer][1]) + w_i[layer][0]) == 0.0)
    {
    };
#else
    while ((res = gsl_ran_gaussian(gsl_gen, w_i[layer][1]) + w_i[layer][0]) == 0.0)
    {
    };
#endif
    return res;
}

double generate_initial_potential(MicrocircuitLayer layer)
{
    double res;
#ifdef MULTIPROCESSING
    while ((res = gsl_ran_gaussian(gsl_gen[0], u_init[layer][1]) + u_init[layer][0]) >= 0.0)
    {
    };
#else
    while ((res = gsl_ran_gaussian(gsl_gen, u_init[layer][1]) + u_init[layer][0]) >= 0.0)
    {
    };
#endif
    return res;
}

uint32_t get_pseudorandom_int(uint32_t start, uint32_t stop)
{
#ifdef MULTIPROCESSING
    return (uint32_t)gsl_ran_flat(gsl_gen[0], (double)start, (double)stop);
#else
    return (uint32_t)gsl_ran_flat(gsl_gen, (double)start, (double)stop);
#endif
}

uint32_t generate_thalamic_spikes(uint8_t thread_id, MicrocircuitLayer layer)
{
#ifdef MULTIPROCESSING
    return gsl_ran_poisson(gsl_gen[thread_id], TIMESTEP * F_TH * thalamic_sizes[layer]);
#else
    return gsl_ran_poisson(gsl_gen, TIMESTEP * F_TH * thalamic_sizes[layer]);
#endif
}

void random_test()
{
#ifdef MULTIPROCESSING
    gsl_gen[0] = gsl_rng_alloc(gsl_rng_mt19937); // allocating and init MT PRNG (seed 0)
    gsl_rng_set(gsl_gen[0], 0);
#else
    gsl_gen = gsl_rng_alloc(gsl_rng_mt19937); // allocating and init MT PRNG (seed 0)
    gsl_rng_set(gsl_gen, 0);
#endif
    uint32_t test_size = 100000;
    double test[test_size];
    double avg = 0.0;
    double variance = 0.0;
    double std = 0.0;

    for (uint32_t i = 0; i < test_size; ++i)
    {
        test[i] = get_pseudorandom_int(0, pop_sizes[0]);
        avg += test[i];
    }
    avg /= test_size;
    for (uint32_t i = 0; i < test_size; ++i)
    {
        variance += pow(test[i] - avg, 2);
    }
    variance /= test_size;
    std = sqrt(variance);
    printf("FLAT AVG: %lf STD: %lf\n", avg, std);
    fflush(stdout);

    variance = 0.0;
    std = 0.0;

    for (uint8_t j = 0; j < 8; ++j)
    {
        avg = 0.0;
        for (uint32_t i = 0; i < test_size; ++i)
        {
            avg += (double)generate_thalamic_spikes(0, j);
        }
        printf("POISSON LAY%u: %lf\n", j, avg / test_size);
        fflush(stdout);
    }
    avg = 0.0;
    variance = 0.0;
    std = 0.0;

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
        num += max_synapses_per_layer[i][layer];
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

void liffifo_init(LIFFifo *fifo)
{
    fifo->num_elements = 0;
    fifo->spikes[0] = UINT32_MAX;
    fifo->spikes[1] = UINT32_MAX;
}

void liffifo_push(LIFFifo *fifo, uint32_t val)
{
    fifo->spikes[fifo->num_elements] = val;
    fifo->num_elements++;
}

uint32_t liffifo_pop(LIFFifo *fifo)
{
    uint32_t holder;
    holder = fifo->spikes[0];
    if (fifo->num_elements)
    {
        switch (fifo->num_elements)
        {
        case 2: // shift, zero 1
            fifo->spikes[0] = fifo->spikes[1];
            fifo->spikes[1] = UINT32_MAX;
            break;
        case 1: // no shift, zero 0
            fifo->spikes[0] = UINT32_MAX;
            break;
        }
        fifo->num_elements--;
    }
    return holder;
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
    liffifo_init(&(neuron->spike_timestamps));
    neuron->total_current = 0.0;
    neuron->presynaptic_current = 0.0;
    neuron->spike_number = 0;
    neuron->thalamic_spikes = 0.0;
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
            if (!max_synapses_per_layer[post_layer][pre_layer])
            {
                printf("FAILED, no synapses between those populations. Skipping.\n");
                continue;
            }
            for (uint32_t sample = 0; sample < max_synapses_per_layer[post_layer][pre_layer]; ++sample)
            {
                synapses[synapse_index].pre_index = get_pseudorandom_int(0, pop_sizes[pre_layer]);
                synapses[synapse_index].post_index = get_pseudorandom_int(0, pop_sizes[post_layer]);
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

void update_neuron(LIFNeuron *neuron, uint32_t current_timestep)
{

    // update membrane
    if (neuron->refractory == TAU_REF)
    {
        liffifo_push(&(neuron->spike_timestamps), current_timestep + neuron->delay);
        neuron->membrane = U_REST;
        neuron->spike = 1;
        neuron->spike_number++;
    }
    else if (neuron->refractory > 0)
    {
        neuron->spike = 0;
    }
    else
        neuron->membrane = P_22 * neuron->membrane + P_21 * neuron->total_current; // P_21 [s/F]

    neuron->total_current *= P_11;                                  // decay the current from previous iteration
    neuron->total_current += neuron->thalamic_spikes * W_EXT * W_F; // add thalamic (input)

    neuron->total_current += neuron->presynaptic_current * W_F;
    // sum = spike_ij * weight_ij [V] [F/s] * [V] = [A]
    // update the total current with thalamics
    // decay rate [1] - this only applies to the "old" spikes. Thalamics are from THIS step
    // thalamic_spikes * W_EXT * W_F; // [F/s] * [V] = [A]
    if (neuron->membrane >= U_THR)
    {
        neuron->refractory = TAU_REF;
    }
    else if (neuron->refractory > 0)
        neuron->refractory--;
    else
        neuron->refractory = 0;
}

void update_network(LIFNetwork *network, uint32_t current_timestep)
{
    uint32_t neuron, pre_neuron;
    LIFNeuron *current_neuron, *pre_neuron_ptr;
    network->current_timestep = current_timestep;

    uint32_t thread_id = 0;

#ifdef MULTIPROCESSING
#pragma omp parallel shared(network) private(thread_id)
    {
        thread_id = omp_get_thread_num();
        gsl_rng_set(gsl_gen[thread_id], thread_id + 123);
#pragma omp barrier
#pragma omp for private(pre_neuron, neuron, pre_neuron_ptr, current_neuron) schedule(runtime)

#endif
        for (neuron = 0; neuron < NEURON_NUMBER; ++neuron)
        {
            current_neuron = &(network->neurons[neuron]);
            current_neuron->presynaptic_current = 0.0; // this needs to be zeroed!

            current_neuron->thalamic_spikes = generate_thalamic_spikes(thread_id, current_neuron->layer);
            for (pre_neuron = 0; pre_neuron < network->neurons[neuron].synapse_count; ++pre_neuron)
            {
                pre_neuron_ptr = current_neuron->presynaptic_neurons[pre_neuron];
                if (network->current_timestep == pre_neuron_ptr->spike_timestamps.spikes[0])
                {
                    current_neuron->presynaptic_current += pre_neuron_ptr->synaptic_amp;
                }
            }
        }

#ifdef MULTIPROCESSING
    }
#endif

#ifdef MULTIPROCESSING
#pragma omp parallel for private(neuron) shared(network) schedule(runtime)
#endif
    for (neuron = 0; neuron < NEURON_NUMBER; ++neuron)
    {
        if (network->neurons[neuron].spike_timestamps.spikes[0] == current_timestep)
        {
            liffifo_pop(&(network->neurons[neuron].spike_timestamps));
        }
        update_neuron(&(network->neurons[neuron]), network->current_timestep);
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

#ifdef MULTIPROCESSING
    for (uint8_t i = 0; i < 16; ++i)
    {
        gsl_gen[i] = gsl_rng_alloc(gsl_rng_mt19937); // allocating and init MT PRNG (seed 0)
    }
    gsl_rng_set(gsl_gen[0], 321); // setting to a known seed for the serial stuff in the beginning
#else
    gsl_gen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(gsl_gen, 321);
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