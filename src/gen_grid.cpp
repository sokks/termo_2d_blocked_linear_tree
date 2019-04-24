#include "grid/grid.h"

int main(int argc, char **argv) {
    int base_level = 0;
    int max_level  = 10;
    string filename = "data/refine/start_grid.dat";
    int n_procs = 1;
    string offsets_filename = "data/refine/start_offsets.dat";
    
    if (argc >= 6) {
        base_level = atoi(argv[1]);
        max_level  = atoi(argv[2]);
        filename   = argv[3];
        n_procs    = atoi(argv[4]);
        offsets_filename   = argv[5];
    }
    GridInit(base_level, max_level);
    BlockedLinearTree grid = BlockedLinearTree(&Area::T0);
    while (grid.MarkToRefine()) {
        grid.DoRefine();
    }
    grid.BuildBlocks();

    grid.Write(filename);
    grid.WriteOffsets(offsets_filename, n_procs);
}