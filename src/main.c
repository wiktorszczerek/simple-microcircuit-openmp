#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "microcircuit.h"
#include "constants.h"

// #define DEBUG

int process_input_args(int argc, char **argv)
{
    if (argc != 2)
    {
        printf("Incorrect number of arguments provided. Expected 1 argument.\n");
        printf("Help: %s -h\n", argv[0]);
        return -1;
    }
    else
    {
        if (!strcmp(argv[1], "-h"))
        {
            puts(help_msg);
            return 0;
        }
        else if (!strcmp(argv[1], "-g"))
            return 1;
        else if (!strcmp(argv[1], "-r"))
            return 2;
    }
    return -2; // don't know what would need to have happened.
}

int main(int argc, char **argv)
{
    int err = process_input_args(argc, argv);
    if (err < 1)
    {
        exit(0);
    }

    LIFNetwork network;

    if (err == 1)
    {
        printf("Not implemented yet!\n");
        save_network(&network); // for now it does nothing!
    }
    else
    {
        // For now this step is run always.
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
    }
    return 0;
}
