#ifndef __MICROCIRCUIT_H
#define __MICROCIRCUIT_H

#include <inttypes.h>

typedef struct LIFNeuron LIFNeuron;
typedef struct LIFConnection LIFConnection;
typedef struct LIFNetwork LIFNetwork;

typedef enum
{
    MLAYER_L23EXC,
    MLAYER_L23INH,
    MLAYER_L4EXC,
    MLAYER_L4INH,
    MLAYER_L5EXC,
    MLAYER_L5INH,
    MLAYER_L6EXC,
    MLAYER_L6INH
} MicrocircuitLayer;

struct LIFNeuron
{
    uint8_t layer;
    float membrane;
    float synaptic_amp;
    uint32_t delay;
    uint8_t refractory;
    uint8_t spike;
    uint32_t spike_number;
    uint32_t spike_timestamps[2];
    uint8_t spike_timestamp_flag;
    uint32_t synapse_count;
    float presynaptic_current;
    float total_current;
    LIFNeuron **presynaptic_neurons;
};

struct LIFConnection
{
    uint8_t pre_layer;
    uint32_t pre_index;
    uint8_t post_layer;
    uint32_t post_index;
};

struct LIFNetwork
{
    LIFNeuron *neurons;
    uint32_t synapse_count;
    uint32_t current_timestep;
};

uint32_t synaptic_number_check();

void initialize_network(LIFNetwork *network);
void deinitialize_network(LIFNetwork *network);

void update_network(LIFNetwork *network, uint32_t current_timestep);
void save_spikes(LIFNetwork *network, uint32_t step);
void save_spiking_rates(LIFNetwork *network);

#endif
