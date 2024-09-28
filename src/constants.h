
/*
    general parameters
*/
// #include <math.h>

// #define DEBUG
#define MULTIPROCESSING

static const float SIMULATION_TIME = 0.01; // s
static const float TIMESTEP = 0.0001;      // s

#ifdef DEBUG
static const int SIMULATION_STEPS = 2;
#else
static const int SIMULATION_STEPS = (int)((float)SIMULATION_TIME / TIMESTEP);
#endif

static const float CAPACITANCE = 0.00000025; // pF
static const float TAU_M = 10;
static const float TAU_REF = 2;
static const float TAU_SYN = 0.5;
static const float U_REST = -65.0;
static const float U_THR = -50.0;
static const float F_TH = 8.0;
static const float W_EXT = 0.15;
static const float W_F = 585.0;    // synapse amp scale
static const float P_11 = 0.82;    // membrane decay
static const float P_22 = 0.99;    // presynaptic decay
static const float P_21 = 0.00036; // injection scale

static const int NEURON_NUMBER = 77169;
static const int LAYER_NUMBER = 8;
static const int INITIAL_SYNAPSE_NUMBER = 20000; // this IS too much and it's okay

#define L23_EXC_POP_SIZE 20683
#define L23_INH_POP_SIZE 5834
#define L4_EXC_POP_SIZE 21915
#define L4_INH_POP_SIZE 5479
#define L5_EXC_POP_SIZE 4850
#define L5_INH_POP_SIZE 1065
#define L6_EXC_POP_SIZE 14395
#define L6_INH_POP_SIZE 2948

//  THALAMIC
#define L23_EXC_THALAMIC_NUMBER 1600
#define L23_INH_THALAMIC_NUMBER 1500
#define L4_EXC_THALAMIC_NUMBER 2100
#define L4_INH_THALAMIC_NUMBER 1900
#define L5_EXC_THALAMIC_NUMBER 2000
#define L5_INH_THALAMIC_NUMBER 1900
#define L6_EXC_THALAMIC_NUMBER 2900
#define L6_INH_THALAMIC_NUMBER 2100

static const int thalamic_sizes[8] =
    {
        L23_EXC_THALAMIC_NUMBER,
        L23_INH_THALAMIC_NUMBER,
        L4_EXC_THALAMIC_NUMBER,
        L4_INH_THALAMIC_NUMBER,
        L5_EXC_THALAMIC_NUMBER,
        L5_INH_THALAMIC_NUMBER,
        L6_EXC_THALAMIC_NUMBER,
        L6_INH_THALAMIC_NUMBER};

static const int pop_sizes[8] =
    {
        L23_EXC_POP_SIZE,
        L23_INH_POP_SIZE,
        L4_EXC_POP_SIZE,
        L4_INH_POP_SIZE,
        L5_EXC_POP_SIZE,
        L5_INH_POP_SIZE,
        L6_EXC_POP_SIZE,
        L6_INH_POP_SIZE};

static const float u_init[8][2] =
    {
        {-68.28, 5.36},
        {-63.16, 4.57},
        {-63.33, 4.74},
        {-63.45, 4.94},
        {-63.11, 4.94},
        {-61.66, 4.55},
        {-66.72, 5.46},
        {-61.43, 4.48}};

static const float w_i[8][2] =
    {
        {0.15, 0.015},
        {-0.6, 0.06},
        {0.15, 0.015},
        {-0.6, 0.06},
        {0.15, 0.015},
        {-0.6, 0.06},
        {0.15, 0.015},
        {-0.6, 0.06}};

static const float delta_i[8][2] =
    {
        {1.5, 0.75},
        {0.75, 0.325},
        {1.5, 0.75},
        {0.75, 0.325},
        {1.5, 0.75},
        {0.75, 0.325},
        {1.5, 0.75},
        {0.75, 0.325}};

static const double synaptic_probability[8][8] =
    {
        {0.1009, 0.1689, 0.0437, 0.0818, 0.0323, 0.0, 0.0076, 0.0},
        {0.1346, 0.1371, 0.0316, 0.0515, 0.0755, 0.0, 0.0042, 0.0},
        {0.0077, 0.0059, 0.0497, 0.1350, 0.0067, 0.0003, 0.0453, 0.0},
        {0.0691, 0.0029, 0.0794, 0.1597, 0.0033, 0.0, 0.1507, 0.0},
        {0.1004, 0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0.0},
        {0.0548, 0.0269, 0.0257, 0.0022, 0.0600, 0.3158, 0.0086, 0.0},
        {0.0156, 0.0066, 0.0211, 0.0166, 0.0572, 0.0197, 0.0396, 0.2252},
        {0.0364, 0.001, 0.0034, 0.0005, 0.0277, 0.008, 0.0658, 0.1443}};

static const int max_synapses_per_layer[8][8] =
    {
        {43163657, 20380255, 19807810, 9269753, 3240096, 0, 2262762, 0},
        {16241459, 4666275, 4040127, 1646172, 2136265, 0, 352718, 0},
        {3490164, 754328, 23869282, 16209759, 712128, 7002, 14290630, 0},
        {7830562, 92698, 9533740, 4794105, 87692, 0, 11885740, 0},
        {10071381, 1759943, 5367532, 151467, 1954720, 1924573, 1424242, 0},
        {1207102, 167136, 599825, 12838, 309915, 358189, 131844, 0},
        {4644616, 554271, 6656342, 1309246, 3993461, 302015, 8205755, 9556691},
        {2219435, 17199, 219659, 8077, 396050, 25117, 2792320, 1254069}};

/*
    Messages.
*/
static const char *help_msg =
    "Possible options are:\n\t-h: View this help.\n\t-g: Generate new microcircuit and store it in the file.\n\t-r: Run the simulation";