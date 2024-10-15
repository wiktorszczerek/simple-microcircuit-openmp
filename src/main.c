#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "microcircuit.h"
#include "constants.h"

double gettime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}

int process_input_args(int argc, char **argv)
{
    if (argc > 3)
    {
        printf("Incorrect number of arguments provided. Expected an option and a specifier.\n");
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
    double dtime, alltime = 0.0;

    if (err == 1)
    {
        printf("Not implemented yet!\n");
    }
    else if (argc != 3)
    {
        printf("Provide number of steps to simulate!");
        exit(1);
    }
    {
        uint32_t sim_steps = (uint32_t)strtol(argv[2], NULL, 10);
        initialize_network(&network);
        printf("Initialization done.\n\n");
        fflush(stdout);
        for (int i = 0; i < sim_steps; ++i)
        {
            printf("Simulating step %d of %d...\r", i + 1, sim_steps);
            fflush(stdout);
            dtime = gettime();
            update_network(&network, i);
            dtime = gettime() - dtime;
            alltime += dtime;
            save_spikes(&network, i);
        }
        alltime /= sim_steps;
        // dtime = gettime() - dtime;
        printf("\nTime of the most recent loop:\t%9.5f seconds\nAverage after %d steps:\t%9.5f seconds\n", dtime, sim_steps, alltime);
        save_spiking_rates(&network);
        deinitialize_network(&network);
    }
    return 0;
}
