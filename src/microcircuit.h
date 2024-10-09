#ifndef __MICROCIRCUIT_H
#define __MICROCIRCUIT_H

#include <inttypes.h>

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

typedef struct LIFNeuronLocation
{
    uint8_t layer;
    uint32_t index;
} LIFNeuronLocation;

typedef struct LIFNeuron
{
    uint8_t layer;
    float membrane;
    float synaptic_amp;
    float delay;
    float refractory;
    uint8_t spike;
    uint32_t synapse_count;
    float presynaptic_current;
    LIFNeuronLocation *presynaptic_neurons;
} LIFNeuron;

typedef struct LIFConnection
{
    uint8_t pre_layer;
    uint32_t pre_index;
    uint8_t post_layer;
    uint32_t post_index;
} LIFConnection;

typedef struct LIFNetwork
{
    LIFNeuron *neurons;
    LIFConnection *synapses;
    uint32_t synapse_count;
    // LIFConnection ***synapses; // 3D? -> 8 x 8 array of arrays
    // LIFNeuron **layers;        // 2D -> 8 layers by howevery many neurons
} LIFNetwork;

uint32_t synaptic_number_check();

void initialize_network(LIFNetwork *network);
void deinitialize_network(LIFNetwork *network);

void update_network(LIFNetwork *network);
void save_spikes(LIFNetwork *network, uint32_t step);

#endif
