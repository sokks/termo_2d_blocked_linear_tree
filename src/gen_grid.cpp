#include "grid/grid.h"

int main(int argc, char **argv) {
    int base_level = 0;
    int max_level  = 10;
    int base_blk_level = 2;
    int max_blk_level  = 5;
    string filename = "data/refine/start_grid.dat";
    int n_procs = 1;
    string blocks_filename = "data/refine/base_grid_blocks.dat";
    string offsets_filename = "data/refine/start_offsets.dat";

    
    if (argc >= 9) {
        base_level = atoi(argv[1]);
        max_level  = atoi(argv[2]);
        base_blk_level = atoi(argv[3]);
        max_blk_level  = atoi(argv[4]);
        filename   = argv[5];
        n_procs    = atoi(argv[6]);
        blocks_filename   = argv[7];
        offsets_filename   = argv[8];
    }
    GridInit(base_level, max_level, base_blk_level, max_blk_level);

    std::cout << "hello 1\n";
    BlockedLinearTree grid(&Area::T0);
    while (grid.MarkToRefine()) {
        std::cout << "hello 2\n";
        grid.RefineBlocks();
        grid.RefineCells();
    }
    std::cout << "hello\n";
    grid.BuildNeighs();

    grid.Decompose(n_procs);

    grid.Write(filename);

    grid.WriteBlocks(blocks_filename);
    grid.WriteOffsets(offsets_filename);
}