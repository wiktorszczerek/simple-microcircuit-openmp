#include "microcircuit.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L23_EXC_POP_SIZE 20683
#define L23_INH_POP_SIZE 5834
#define L4_EXC_POP_SIZE 21915
#define L4_INH_POP_SIZE 5479
#define L5_EXC_POP_SIZE 4850
#define L5_INH_POP_SIZE 1065
#define L6_EXC_POP_SIZE 14395
#define L6_INH_POP_SIZE 2948

#define NEURON_NUMBER   77169

const int pop_sizes[8] =
{
  L23_EXC_POP_SIZE,
  L23_INH_POP_SIZE,
  L4_EXC_POP_SIZE,
  L4_INH_POP_SIZE,
  L5_EXC_POP_SIZE,
  L5_INH_POP_SIZE,
  L6_EXC_POP_SIZE,
  L6_INH_POP_SIZE
};


const float u_init[8][2] =
{
  {-68.28, 5.36},
  {-63.16, 4.57},
  {-63.33, 4.74},
  {-63.45, 4.94},
  {-63.11, 4.94},
  {-61.66, 4.55},
  {-66.72, 5.46},
  {-61.43, 4.48}
};

const float w_i[8][2] =
{
  {0.15, 0.015},
  {-0.6, 0.06},
  {0.15, 0.015},
  {-0.6, 0.06},
  {0.15, 0.015},
  {-0.6, 0.06},
  {0.15, 0.015},
  {-0.6, 0.06}
};

const float delta_i[8][2] =
{
  {1.5, 0.75},
  {0.75, 0.325},
  {1.5, 0.75},
  {0.75, 0.325},
  {1.5, 0.75},
  {0.75, 0.325},
  {1.5, 0.75},
  {0.75, 0.325}
};

const float synaptic_probability[8][8] = 
{
  {0.1009, 0.1689, 0.0437, 0.0818, 0.0323, 0.0, 0.0076, 0.0},
  {0.1346, 0.1371, 0.0316, 0.0515, 0.0755, 0.0, 0.0042, 0.0},
  {0.0077, 0.0059, 0.0497, 0.1350, 0.0067, 0.0003, 0.0453, 0.0},
  {0.0691, 0.0029, 0.0794, 0.1597, 0.0033, 0.0, 0.1507, 0.0},
  {0.1004, 0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0.0},
  {0.0548, 0.0269, 0.0257, 0.0022, 0.0600, 0.3158, 0.0086, 0.0},
  {0.0364, 0.001, 0.0034, 0.0005, 0.0277, 0.008, 0.0658, 0.1443}
};

// TODO: add a license here
double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
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

void create_neuron(LIFNeuron* neuron, MicrocircuitLayer layer)
{
  neuron->membrane = generate_initial_potential(layer);
  neuron->synaptic_amp = generate_synaptic_amp(layer);
  neuron->delay = generate_delay(layer);
}

void generate_synapse(LIFNetwork* network, MicrocircuitLayer pre_layer, MicrocircuitLayer post_layer, int pre_index, int post_index)
{
  network->synapses = (LIFSynapse*)realloc(network->synapses, sizeof(network->synapses) + sizeof(LIFSynapse));
  network->synapses->pre_layer = pre_layer;
  network->synapses->post_layer = post_layer;
  network->synapses->pre_index = pre_index;
  network->synapses->post_index = post_index;
}

float generate_uniform_probability()
{
  return (float)rand()/RAND_MAX;
}

void initialize_synapses(LIFNetwork* network)
{
    network->synapses = (LIFSynapse*)malloc(sizeof(LIFSynapse));
}

void create_synapses(LIFNetwork* network, MicrocircuitLayer pre_layer)
{
  // int counter = 0;
  //TODO: those loops can be optimized (indexing)
  for(int pre_neuron=0; pre_neuron < pop_sizes[pre_layer]; ++pre_neuron)
  {
    for(int post_layer=0; post_layer < 8; ++post_layer)
    {
      for(int post_neuron=0; post_neuron < pop_sizes[post_layer]; ++post_neuron)
      {
        if(generate_uniform_probability() <= synaptic_probability[pre_layer][post_layer])
        {
          // if(pre_layer == MLAYER_L23INH && post_layer == MLAYER_L5EXC)
          //   counter++;
          generate_synapse(network, pre_layer, post_layer, pre_neuron, post_neuron);
        }
      }
    }
  }
  // printf("Counter L23INH - L5EXC = %d\n", counter);

  // for(int i=0; i<8; ++i)
  // {
  //   for(int j=0; j<pop_sizes[i]; ++j)
  //   {
  //     if(generate_probability() <= synaptic_probability[layer][i])
  //     {
  //       if(network->synapses == NULL)
  //         network->synapses = (LIFSynapse*)malloc(sizeof(LIFSynapse));
  //       else
  //         network->synapses = (LIFSynapse*)realloc(network->synapses, sizeof(network->synapses) + sizeof(LIFSynapse));
  //       generate_synapse(network, layer, i, j, )
  //     }
  //   }
  // }
}

void initialize_layer(LIFNetwork* network, MicrocircuitLayer layer)
{
  int pop_size = 0;
  LIFNeuron* layer_neuron;
  switch(layer)
  {
    case MLAYER_L23EXC:
      network->l23_exc = (LIFNeuron*)malloc(L23_EXC_POP_SIZE*sizeof(LIFNeuron));
      pop_size = L23_EXC_POP_SIZE;
      layer_neuron = network->l23_exc;
      break;
    case MLAYER_L23INH:
      network->l23_inh = (LIFNeuron*)malloc(L23_INH_POP_SIZE*sizeof(LIFNeuron));
      pop_size = L23_INH_POP_SIZE;
      layer_neuron = network->l23_inh;
      break;
    case MLAYER_L4EXC:
      network->l4_exc = (LIFNeuron*)malloc(L4_EXC_POP_SIZE*sizeof(LIFNeuron));
      pop_size = L4_EXC_POP_SIZE;
      layer_neuron = network->l4_exc;
      break;
    case MLAYER_L4INH:
      network->l4_inh = (LIFNeuron*)malloc(L4_INH_POP_SIZE*sizeof(LIFNeuron));
      pop_size = L4_INH_POP_SIZE;
      layer_neuron = network->l4_inh;
      break;
    case MLAYER_L5EXC:
      network->l5_exc = (LIFNeuron*)malloc(L5_EXC_POP_SIZE*sizeof(LIFNeuron));
      pop_size = L5_EXC_POP_SIZE;
      layer_neuron = network->l5_exc;
      break;
    case MLAYER_L5INH:
      network->l5_inh = (LIFNeuron*)malloc(L5_INH_POP_SIZE*sizeof(LIFNeuron));
      pop_size = L5_INH_POP_SIZE;
      layer_neuron = network->l5_inh;
      break;
    case MLAYER_L6EXC:
      network->l6_exc = (LIFNeuron*)malloc(L6_EXC_POP_SIZE*sizeof(LIFNeuron));
      pop_size = L6_EXC_POP_SIZE;
      layer_neuron = network->l6_exc;
      break;
    case MLAYER_L6INH:
      network->l6_inh = (LIFNeuron*)malloc(L6_INH_POP_SIZE*sizeof(LIFNeuron));
      pop_size = L6_INH_POP_SIZE;
      layer_neuron = network->l6_inh;
      break;
    default:
      printf("Can't initialize layer - layer value unknown: %d\n", layer);
      exit(1);
  }

  for(int i=0; i<pop_size; ++i)
    create_neuron(layer_neuron, layer);

  initialize_synapses(network);
  create_synapses(network, layer);

}



void initialize_network(LIFNetwork* network) {

  // in here, for the parallel code the workeer will realize one layer in generation
  for(int i=0; i<8; ++i)
  {
    printf("Initializing Layer %d...", i);
    fflush(stdout);
    initialize_layer(network, i); 
    printf("Done!\n");
    fflush(stdout);
  }

  // neurons = (LIFNeuron*)malloc(NEURON_NUMBER * sizeof(LIFNeuron));



  // network->l23_exc =
  //     (LIFNeuron *)malloc(L23_EXC_POP_SIZE * sizeof(LIFNeuron));


  // for (int i = 0; i < L23_EXC_POP_SIZE; ++i) {
  //   network->l23_exc[i].membrane = generate_initial_potential(MLAYER_L23EXC);
  //   network->l23_exc[i].synaptic_amp = generate_synaptic_amp(MLAYER_L23EXC);
  //   network->l23_exc[i].delay = generate_delay(MLAYER_L23EXC);
  // }

  // for(int i = 0; i < L23_EXC_POP_SIZE; ++ i)
  // {
  //   for(int j=0; j < 8; ++j)
  //   {

  //   }


  //   if(((float)rand())/RAND_MAX <= synaptic_probability[i][j]
  //   network->l23_exc[i]
  // }
}

void deinitialize_network(LIFNetwork* network)
{
  free(network->l23_exc);
  free(network->l23_inh);
  free(network->l4_exc);
  free(network->l4_inh);
  free(network->l5_exc);
  free(network->l5_inh);
  free(network->l6_exc);
  free(network->l6_inh);
  free(network->synapses);
}