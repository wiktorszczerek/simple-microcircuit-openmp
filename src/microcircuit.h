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

typedef struct LIFSynapse
{
    int pre_layer;
    int post_layer;
    int pre_index;
    int post_index;

} LIFSynapse;

typedef struct LIFNeuron
{
    float membrane;
    float synaptic_amp;
    float delay;

} LIFNeuron;

typedef struct LIFNetwork
{
    LIFNeuron* l23_exc;
    LIFNeuron* l23_inh;
    LIFNeuron* l4_exc;
    LIFNeuron* l4_inh;
    LIFNeuron* l5_exc;
    LIFNeuron* l5_inh;
    LIFNeuron* l6_exc;
    LIFNeuron* l6_inh;
    LIFSynapse* synapses;
} LIFNetwork;

void initialize_network(LIFNetwork* network);
void deinitialize_network(LIFNetwork* network);


#endif
