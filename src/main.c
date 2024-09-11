#include <stdio.h>
#include <stdlib.h>
#include "microcircuit.h"


int main(int argc, char **argv) { 
    
    LIFNetwork network;
    initialize_network(&network);

    deinitialize_network(&network);
    return 0; 
}
