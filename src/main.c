#include <stdio.h>
#include <stdlib.h>
#include "microcircuit.h"
#include "constants.h"

int main(int argc, char **argv)
{
    LIFNetwork network;
    initialize_network(&network);
    printf("Initialization done.\n\n");
    fflush(stdout);
    for (int i = 0; i < SIMULATION_STEPS; ++i)
    {
        printf("Simulating step %d of %d...\r", i + 1, SIMULATION_STEPS);
        fflush(stdout);
        update_network(&network);
        save_spikes(&network, i);
    }

    deinitialize_network(&network);
    return 0;
}
