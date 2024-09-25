#ifndef __MICROCIRCUIT_H
#define __MICROCIRCUIT_H

#include <inttypes.h>
// #include "constants.h"

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

// typedef struct LIFSynapse
// {
//     int pre_layer;
//     // int post_layer;
//     int pre_index;
//     // int post_index;
// } LIFSynapse;

typedef struct LIFNeuron
{
    float membrane;
    float synaptic_amp;
    float delay;
    float refractory;
    uint8_t spike;
    // int (*synapses)[2]; // postsynaptic connections
    int synapse_count;
    float presynaptic_current;
} LIFNeuron;

typedef struct LIFConnection
{
    int pre_index;
    int post_index;
} LIFConnection;

typedef struct LIFNetwork
{
    LIFConnection ***synapses; // 3D? -> 8 x 8 array of arrays
    LIFNeuron **layers;        // 2D -> 8 layers by howevery many neurons
    uint8_t *activity;         // 1D
} LIFNetwork;

int synaptic_number_check();

void initialize_network(LIFNetwork *network);
void deinitialize_network(LIFNetwork *network);

void update_network(LIFNetwork *network);
void save_network(LIFNetwork *network);
void save_spikes(LIFNetwork *network, int step);

#endif
