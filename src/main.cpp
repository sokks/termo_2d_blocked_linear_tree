#include "proc/proc.h"

using std::string;
using std::cout;
using std::endl;


bool WRITE_LAYERS = true;
int write_freq = 100;
string baseFolderTemp = "data/temp/";

string gen_filename(string baseFolder, int n) {
    string num = std::to_string(n);
    int max = 6;
    int additional = max - num.length();
    string res = string(additional, '0') + num;
    string fname = baseFolder + res + ".out";
    return fname;
}


int main(int argc, char **argv) {
    if (argc < 6) {
        std::cout << "usage: prog <base_lvl> <max_lvl> <offsets_file> <grid_file> <time_steps>\n";
        return 0;
    }
    int    base_level     = atoi(argv[1]);
    int    max_level      = atoi(argv[2]);
    int    base_blk_level = atoi(argv[3]);
    int    max_blk_level  = atoi(argv[4]);
    string offsets_file   = argv[5];
    string grid_file      = argv[6];
    int    ts_n           = std::atoi(argv[7]);

    GridInit(base_level, max_level, base_blk_level, max_blk_level);
    SolverInit(ts_n);

    Proc p;
    p.MPIInit(argc, argv);

    cout << "base_level=" << base_level << 
            " max_level=" << max_level <<
            " offsets_file=" << offsets_file << 
            " grid_file=" << grid_file << 
            " ts_n=" << ts_n << endl;
    
    p.InitMesh(offsets_file, grid_file, &Area::T0);

    // usleep(10000000);

     p.BuildGhosts();

    // usleep(10000000);

     for (int k = 0; k < ts_n; k++) {
         if (WRITE_LAYERS && (k%write_freq ==0)) {
             p.WriteT(gen_filename(baseFolderTemp, k));
         }
         p.MakeStep();
     }

    p.WriteT(gen_filename(baseFolderTemp, ts_n));
    p.WriteStat("data/stat.out");

    p.MPIFinalize();
}