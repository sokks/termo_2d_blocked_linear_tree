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
        base_blk_lvl = atoi(argv[3]);
        max_blk_lvl  = atoi(argv[4]);
        filename   = argv[5];
        n_procs    = atoi(argv[6]);
        offsets_filename   = argv[7];
    }
    GridInit(base_level, max_level, base_blk_lvl, max_blk_lvl);

    BlockedLinearTree grid(&Area::T0);
    // int a = 1;
    while (grid.MarkToRefine()) {
        grid.RefineBlocks();
        grid.RefineCells();
        // a++;
    }

    grid.Decompose(n_procs);

    grid.Write(filename);
    // grid.WriteOffsets(offsets_filename, n_procs);
}