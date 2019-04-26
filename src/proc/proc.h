#include <vector>
#include <mpi.h>
#include <math.h>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <unistd.h>

#include "../area/area.h"
#include "../grid/grid.h"

using std::vector;
using std::map;
using std::string;

extern double base_tau;
extern double tau;
extern int    time_steps;

void SolverInit(int ts_n);


/* * * * * * * * * * * СТАТИСТИКА * * * * * * * * * * */

class MpiTimer {
    double s_time;
    double dur = 0;

public:
    void   Start() { s_time = MPI_Wtime(); }
    double Stop() { double cur_dur = MPI_Wtime() - s_time; dur += cur_dur; return cur_dur; }
    double FullDur() { return dur; }
};

struct MpiInfo {
    MPI_Comm comm;
    int comm_rank;
    int comm_size;

    string toString();
};

struct Stat {
    int proc_cells;
    map<string, MpiTimer> timers;

    string toString();
};



/* * * * * * * * * * * РАБОТА ПРОЦЕССА * * * * * * * * * * */

struct MetaInfo {
    GlobalNumber_t procStart;
    GlobalNumber_t procEnd;
    GlobalNumber_t procG;

    string toString();
};

struct FullMeta {
    vector<GlobalNumber_t> metas;
    int one_meta_len = sizeof(GlobalNumber_t) * 3;

    FullMeta(){}

    FullMeta(MetaInfo my_meta, int n_procs, int my_rank) {
        metas = vector<GlobalNumber_t>(n_procs * 3);
        metas[my_rank*3]   = my_meta.procStart;
        metas[my_rank*3+1] = my_meta.procEnd;
        metas[my_rank*3+2] = my_meta.procG;
    }

    MetaInfo GetMetaOfProc(int proc_n) {
        MetaInfo m;
        m.procStart = metas[proc_n*3];
        m.procEnd   = metas[proc_n*3+1];
        m.procG     = metas[proc_n*3+2];
        return m;
    }
};

class Proc {
    MpiInfo mpiInfo;
    Stat    stat;

    int time_step_n = 0;

    
    GlobalNumber_t totalG;
    GlobalNumber_t procG;
    GlobalNumber_t offsetG;
    FullMeta meta;
    BlockedLinearTree mesh;

    vector<GlobalNumber_t> *ghosts_out_ids = nullptr;
    LinearTree *ghosts_in  = nullptr;
    vector<double> *ghosts_in_temps  = nullptr;
    vector<double> *ghosts_out_temps = nullptr;

    int active_neighs_num = 0;

public:
    Proc();
    ~Proc();

    /// MPIInit must be called before start work!
    int MPIInit(int argc, char **argv);
    int MPIFinalize();

    /// InitMesh создает базовую структуру сетки
//    int InitMesh();
    int InitMesh(string offsets_filename, string blocks_filename, double (*start_func)(double, double));

    // InitNeighs прописывает соседей блоков
    int InitNeighs();

    // BuildGhosts строит структуры для обмена границами
    int BuildGhosts();

    // ExchangeGhosts обменивается с соседями
    int StartExchangeGhosts();
    int StopExchangeGhosts();

    // MakeStep делает один временной шаг
    void MakeStep();

    // I/O
    void WriteT(string filename);
    void WriteStat(string filename);

private:

    int find_owner(GlobalNumber_t cell_id);
    vector<char> get_write_cells_struct();

    int ISendGhosts();
    int IRecvGhosts();
    int WaitallGhosts();

    void PrintMyCells();
    void PrintGhostCells();

// запрещаем копирование
    Proc(Proc&);
    Proc& operator=(Proc&);
};

