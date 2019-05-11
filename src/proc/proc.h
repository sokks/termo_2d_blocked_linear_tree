#include <vector>
#include <mpi.h>
#include <math.h>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <unistd.h>
// #include <omp.h>

// #include "../area/area.h"
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
    double dur;

public:
    MpiTimer() { dur = 0; }
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
    int one_meta_len;

    FullMeta(){one_meta_len = sizeof(GlobalNumber_t) * 3;}

    FullMeta(MetaInfo my_meta, int n_procs, int my_rank) {
        one_meta_len = sizeof(GlobalNumber_t) * 3;
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

    FullMeta meta;
    BlockedLinearTree mesh;


    map<GlobalNumber_t, BlockOfCells*> fake_ghost_blocks;  // блоки без ячеек (map[blk_idx]blk)
    vector<vector<GlobalNumber_t> > fake_blocks_out_ids;   // айдишники исходящих блоков (первая размерность - ранг процессора)
    vector<vector<GlobalNumber_t> > fake_blocks_in_ids;    // айдишники входящих блоков  (первая размерность - ранг процессора)

    vector<map<GlobalNumber_t, vector<GlobalNumber_t > > > blocks_cells_in_idxs;  // []map[blk_idx][]blk_cells (первая размерность - ранг процессора)
    vector<map<GlobalNumber_t, vector<GlobalNumber_t > > > blocks_cells_out_idxs;

    vector<vector<int>> blocks_cells_out_lens;       // длины массивов айдишников исходящих ячеек по процессам
    vector<vector<int>> blocks_cells_in_lens;        // длины массивов айдишников входящих блоков по процессам
    vector<vector<GlobalNumber_t > > cells_out_idxs; // айдишники исходящих ячеек по процессам
    vector<vector<GlobalNumber_t > > cells_in_idxs;  // айдишники входящих ячеек по процессам
    vector<vector<double> > cells_out_idxs_temps;    // соответствующие исходящие температуры
    vector<vector<double> > cells_in_idxs_temps;     // соответствующие входящие температуры

    /// active_neights_num is a number of processors that we have dependencies with
    int active_neighs_num = 0;

    vector<MPI_Request> send_reqs;
    vector<MPI_Status>  send_statuses;
    vector<MPI_Request> recv_reqs;
    vector<MPI_Status>  recv_statuses;

public:
    Proc();
    ~Proc();

    /// MPIInit must be called before start work!
    int MPIInit(int argc, char **argv);
    int MPIFinalize();

    /// InitMesh параллельно считывает сетку из файла и инициализирует ячейки
    int InitMesh(string offsets_filename, string blocks_filename, double (*start_func)(double, double));

    /// BuildGhosts строит структуры для обмена границами
    int BuildGhosts();

    /// ExchangeGhosts обменивается с соседями
    int StartExchangeGhosts();
    int StopExchangeGhosts();

    /// MakeStep делает один временной шаг
    void MakeStep();

    /// I/O
    void WriteT(string filename);
    void WriteStat(string filename);

    size_t GetProcAllocMem();

private:

    int find_owner(GlobalNumber_t cell_id);
//    vector<char> get_write_cells_struct();

    void build_fake_ghost_blocks();
    void build_ghost_cells();
    void build_needed_cells_for_blocks_map();

    void get_border_cond(char *cond_type, double (**cond_func)(double, double, double), Neigh border);

    void PrintMyBlocks();
    void PrintGhostCells();

// запрещаем копирование
    Proc(Proc&);
    Proc& operator=(Proc&);
};

