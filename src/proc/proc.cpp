#include "proc.h"

#define FULL_DEBUG 0

#define USE_SIMPLE_FLOWS 1
// #define USE_OPEN_MP 1

using std::vector;
using std::map;
using std::string;
using std::sort;

using std::cout;
using std::endl;

double base_tau = 0.001;
double      tau = 0.001;
int    time_steps = 100;

void SolverInit(int ts_n) {
    base_tau = (base_dx * base_dx) / 4;
    tau = pow((min_dx * min_dx), 2) / 4;
    time_steps = ts_n;

    std::cout << "tau=" << tau << " t_end=" << tau * time_steps << std::endl;
}


/* * * * * * * * * * * СТАТИСТИКА * * * * * * * * * * */

string MpiInfo::toString() {
    std::ostringstream stringStream;
    stringStream << "{ comm: " << comm << ", ";
    stringStream << "comm_size: " << comm_size << ", ";
    stringStream << "comm_rank: " << comm_rank << " }";

    return stringStream.str();
}

string Stat::toString() {
    std::ostringstream stringStream;
    for (map<string, MpiTimer>::iterator it = timers.begin(); it != timers.end(); it++) {
        stringStream << it->first << ": " << it->second.FullDur() << std::endl;
    }
    return stringStream.str();
}

string MetaInfo::toString() {
    std::ostringstream stringStream;
    stringStream << "{ (first: " << procStart << " ";
    stringStream << "last: " << procEnd << "), ";
    stringStream << "len: " << procG << " }";

    return stringStream.str();
}


/* * * * * * * * * * * РАБОТА ПРОЦЕССА * * * * * * * * * * */

Proc::Proc() {
    stat.timers["total"] = MpiTimer();
    stat.timers["io"] = MpiTimer();
    stat.timers["build_ghosts"] = MpiTimer();
    stat.timers["exchange_ghosts"] = MpiTimer();
    stat.timers["communication"] = MpiTimer();
    stat.timers["find_cell"] = MpiTimer();
    stat.timers["step"] = MpiTimer();
    stat.timers["get_border_cond"] = MpiTimer();
    stat.timers["get_possible_neighs"] = MpiTimer();
    stat.timers["sort_neighs"] = MpiTimer();
    stat.timers["compute_temps"] = MpiTimer();
    stat.timers["compute_temps_border"] = MpiTimer();
    stat.timers["compute_temps_border_inblock"] = MpiTimer();
    stat.timers["compute_temps_border_edges"] = MpiTimer();
    stat.timers["compute_temps_border_interblock"] = MpiTimer();
    stat.timers["compute_temps_border_res"] = MpiTimer();

    time_step_n = 0;
    active_neighs_num = 0;
}

Proc::~Proc() {

    cout << mpiInfo.comm_rank << " DESTROYING PROC\n";

    // std::cout << mpiInfo.comm_rank << " ";
    // std::cout << stat.toString() << std::endl;
}

int Proc::MPIInit(int argc, char **argv) {
    // MPI_Init(&argc, &argv);
    int provided;
    #ifdef USE_OPEN_MP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    #else
    MPI_Init(&argc, &argv);
    #endif
    mpiInfo.comm = MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiInfo.comm_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiInfo.comm_size);

    // cout << mpiInfo.comm_rank << " N_THREADS=" << omp_get_num_threads() << endl;

    stat.timers["total"].Start();
    return 0;
}
int Proc::MPIFinalize() {
    stat.timers["total"].Stop();

    MPI_Finalize();
    return 0;
}

int Proc::InitMesh(char* offsets_filename, char* blocks_filename, double (*start_func)(double, double)) {
    stat.timers["init_mesh"] = MpiTimer();
    stat.timers["init_mesh"].Start();

    // [0] - offset, [1] - len
    int range[2];
    int one_sz = 5 * sizeof(int);

    MPI_File fh;
    MPI_File_open( mpiInfo.comm, offsets_filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_at(fh, 2 * mpiInfo.comm_rank * sizeof(int), range, 2, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    std::cout << mpiInfo.comm_rank << " read range[0]=" << range[0] << " range[1]=" << range[1] << std::endl;

    vector<char> buffer(range[1], 1);
    std::cout << "will read cells file\n";

    MPI_File_open( mpiInfo.comm, blocks_filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_at(fh, range[0], &buffer[0], range[1], MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    std::cout << "read cells file buffer.size()=" << buffer.size() << std::endl;
    cout << mpiInfo.comm_rank << " buf[0]=" << int(buffer[0]) << endl;

    mesh.GenFromWriteBlocksStruct(buffer, start_func);

    stat.timers["init_mesh"].Stop();
    std::cout << mpiInfo.comm_rank <<  " MESH INITED   " <<
                 "n_of_blocks=" << mesh.blocks.size() <<
                " blocks_offset=" << range[0] / one_sz << std::endl;
    

    MetaInfo my_meta;
    my_meta.procStart = mesh.blocks[0].idx.get_global_number();
    my_meta.procEnd   = mesh.blocks[mesh.blocks.size()-1].idx.get_global_number();
    my_meta.procG     = mesh.blocks.size();

    cout << mpiInfo.comm_rank << " my_meta={" << my_meta.procStart << "," << my_meta.procEnd << "," << my_meta.procG << "}\n";

    meta = FullMeta(my_meta, mpiInfo.comm_size, mpiInfo.comm_rank);
    // for (int i = 0; i < mpiInfo.comm_size * 3; i++) {
    //     cout << "full_meta[" << i << "]=" << meta.metas[i] << endl;
    // }
    // cout << mpiInfo.comm_rank << " sizeof(my_meta)=" << sizeof(MetaInfo) << "sizeof(metas)=" << meta.metas.size() * meta.one_meta_len << endl;
    MPI_Allgather( (void*)(meta.metas.data() + mpiInfo.comm_rank*3), 3, MPI_LONG_LONG_INT, meta.metas.data(), 3, MPI_LONG_LONG_INT, mpiInfo.comm);

    std::cout << mpiInfo.comm_rank <<  " exchanged meta  ";
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        std::cout << meta.GetMetaOfProc(i).toString() << " | ";
    }
    
    return 0;
}


int Proc::BuildGhosts() {
    MPI_Barrier(mpiInfo.comm);

    cout << mpiInfo.comm_rank << " building ghosts\n";

    stat.timers["build_ghosts"].Start();

    // (1) create fake_ghost_blocks struct
    build_fake_ghost_blocks();

    // (2) create each block req indices
    build_ghost_cells();

    // (3) build border cells neighs pointers
    // TODO
    build_border_cells_neighs_ptrs();

    stat.timers["build_ghosts"].Stop();
    MPI_Barrier(mpiInfo.comm);

    cout << mpiInfo.comm_rank << " built ghosts\n";
    return 0;
}


void Proc::build_fake_ghost_blocks() {

    cout << mpiInfo.comm_rank << " build_fake_ghost_blocks started\n";

    // строим списки индексов соседей блоков по процессам
    fake_blocks_out_ids = vector<vector<GlobalNumber_t> >(mpiInfo.comm_size);
    for (int blk_i = 0; blk_i < mesh.blocks.size(); blk_i++) {
        GlobalNumber_t blk_num = mesh.blocks[blk_i].idx.get_global_number();
        // cout << mpiInfo.comm_rank << " " << blk_num << "  ";
        // if (mesh.blocks[blk_i].idx.get_global_number() == GlobalNumber_t(1536)) {
            // cout << " L: ";
            // for (GlobalNumber_t neight_num: mesh.blocks[blk_i].neighs_left_idxs) {
            //     cout << neight_num << " ";
            // }
            // cout << " R: ";
            // for (GlobalNumber_t neight_num: mesh.blocks[blk_i].neighs_right_idxs) {
            //     cout <<  neight_num << " ";
            // }
            // cout << " D: ";
            // for (GlobalNumber_t neight_num: mesh.blocks[blk_i].neighs_down_idxs) {
            //     cout << neight_num << " ";
            // }
            // cout << " U: ";
            // for (GlobalNumber_t neight_num: mesh.blocks[blk_i].neighs_upper_idxs) {
            //     cout << neight_num << " ";
            // }
            // cout << endl;
        // }
        for (int jjj = 0; jjj < mesh.blocks[blk_i].neighs_left_idxs.size(); jjj++) {
            GlobalNumber_t neight_num = mesh.blocks[blk_i].neighs_left_idxs[jjj];
            int o = find_owner(neight_num);
            if (o != mpiInfo.comm_rank) {
                fake_blocks_out_ids[o].push_back(blk_num);
            }
        }
        for (int jjj = 0; jjj < mesh.blocks[blk_i].neighs_right_idxs.size(); jjj++) {
            GlobalNumber_t neight_num = mesh.blocks[blk_i].neighs_right_idxs[jjj];
            int o = find_owner(neight_num);
            if (o != mpiInfo.comm_rank) {
                fake_blocks_out_ids[o].push_back(blk_num);
            }
        }
        for (int jjj = 0; jjj < mesh.blocks[blk_i].neighs_upper_idxs.size(); jjj++) {
            GlobalNumber_t neight_num = mesh.blocks[blk_i].neighs_upper_idxs[jjj];
            int o = find_owner(neight_num);
            if (o != mpiInfo.comm_rank) {
                fake_blocks_out_ids[o].push_back(blk_num);
            }
        }
        for (int jjj = 0; jjj < mesh.blocks[blk_i].neighs_down_idxs.size(); jjj++) {
            GlobalNumber_t neight_num = mesh.blocks[blk_i].neighs_down_idxs[jjj];
            int o = find_owner(neight_num);
            if (o != mpiInfo.comm_rank) {
                fake_blocks_out_ids[o].push_back(blk_num);
            }
        }
    }

    // cout << "PARAM1\n";

    // сортируем и удаляем повторения
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        std::sort(fake_blocks_out_ids[i].begin(), fake_blocks_out_ids[i].end());
        vector<GlobalNumber_t>::iterator last = std::unique(fake_blocks_out_ids[i].begin(), fake_blocks_out_ids[i].end());
        fake_blocks_out_ids[i].erase(last, fake_blocks_out_ids[i].end());
    }

    // формируем списки длин по процессам для дальнейшего обмена
    vector<int> out_lens(mpiInfo.comm_size, 0);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        out_lens[i] = fake_blocks_out_ids[i].size();
    }

    // здесь определяется количество активных соседей
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        if (out_lens[i] > 0) {
            active_neighs_num++;
        }
    }
    // cout << mpiInfo.comm_rank << " active_neighs_num=" << active_neighs_num << endl;

    send_reqs = vector<MPI_Request>(active_neighs_num);
    send_statuses = vector<MPI_Status>(active_neighs_num);
    recv_reqs = vector<MPI_Request>(active_neighs_num);
    recv_statuses = vector<MPI_Status>(active_neighs_num);

//  cout << "PARAM2\n";
    // обмениваемся длинами чтобы знать сколько блоков от кого принимать
    vector<int> in_lens(mpiInfo.comm_size, 0);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        MPI_Gather(&out_lens[i], 1, MPI_INT, &in_lens[0], 1, MPI_INT, i, mpiInfo.comm);
    }
//  cout << "PARAM3\n";

    // (1) обмен айдишниками блоков
    // формируем массивы для приема
    fake_blocks_in_ids = vector<vector<GlobalNumber_t> >(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        fake_blocks_in_ids[i] = vector<GlobalNumber_t>(in_lens[i]);
    }
//  cout << "PARAM4\n";
    // отправляем и принимаем
    int req_num = 0;
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if (out_lens[n] > 0) {
            MPI_Isend(&fake_blocks_out_ids[n][0], out_lens[n], MPI_LONG_LONG_INT, n, 0, mpiInfo.comm, &send_reqs[req_num]);
            req_num++;
        }
    }
    //  cout << "PARAM5\n";
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ( (n != mpiInfo.comm_rank) && (in_lens[n] > 0) ) {
            MPI_Status status;
            MPI_Recv(&fake_blocks_in_ids[n][0], in_lens[n], MPI_LONG_LONG_INT, n, MPI_ANY_TAG, mpiInfo.comm, &status);
        }
    }

    // cout << mpiInfo.comm_rank << " fake_blocks_in_ids: ";
    // for (int n = 0; n < mpiInfo.comm_size; n++) {
    //     cout << "{" << n << ": ";
    //     for (GlobalNumber_t idx: fake_blocks_in_ids[n]) {
    //         cout << idx << " ";
    //     }
    //     cout << "} ";
    // }
    // cout << endl;
    //  cout << "PARAM6\n";
    MPI_Waitall(active_neighs_num, &send_reqs[0], &send_statuses[0]);


    // (2) обмен информацией о блоках
    // формируем массивы для отправки
    vector<vector<int> > data_fake_blocks_out;
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        data_fake_blocks_out.push_back(vector<int>());
        for (int j = 0; j < fake_blocks_out_ids[i].size(); j++) {
            cout << mpiInfo.comm_rank << " " << fake_blocks_out_ids[i][j] << endl;
            BlockOfCells *blk = mesh.find_block(fake_blocks_out_ids[i][j]);

            data_fake_blocks_out[i].push_back(blk->idx.lvl);
            data_fake_blocks_out[i].push_back(blk->cells_lvl);
            data_fake_blocks_out[i].push_back(blk->sz);
        }
    }
    //  cout << "PARAM7\n";
    // формируем массивы для приема
    vector<vector<int> > data_fake_blocks_in;
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        data_fake_blocks_in.push_back(vector<int>(in_lens[i]*3));
    }

    // cout << "PARAM8\n";
    // отправляем и принимаем
    req_num = 0;
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if (out_lens[n] > 0) {
            MPI_Isend(&data_fake_blocks_out[n][0], out_lens[n] * 3, MPI_INT, n, 0, mpiInfo.comm, &send_reqs[req_num]);
            req_num++;
        }
    }
    // cout << "PARAM9\n";
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ( (n != mpiInfo.comm_rank) && (in_lens[n] > 0) ) {
            MPI_Status status;
            MPI_Recv(&data_fake_blocks_in[n][0], in_lens[n] * 3, MPI_INT, n, MPI_ANY_TAG, mpiInfo.comm, &status);
        }
    }
    MPI_Waitall(active_neighs_num, &send_reqs[0], &send_statuses[0]);

    // cout << "PARAM10\n";
    // (3) формируем собственно мапу
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        for (int j = 0; j < in_lens[i]; j++) {
            fake_ghost_blocks[fake_blocks_in_ids[i][j]] = new BlockOfCells(
                    data_fake_blocks_in[i][j*3+1], // cells_lvl
                    data_fake_blocks_in[i][j*3],   // blk_lvl
                    data_fake_blocks_in[i][j*3+2], // sz
                    fake_blocks_in_ids[i][j]);   // idx
        }
    }
    // cout << "PARAM11\n";

    cout << mpiInfo.comm_rank << " build_fake_ghost_blocks finished\n";
}

void Proc::build_ghost_cells() {

    // (1) построение мапки нужных ячеек по блокам
    build_needed_cells_for_blocks_map();


    // (2) обмен этими нужными ячейками

    // обмен длинами
    blocks_cells_out_lens = vector<vector<int> >(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        for (int j = 0; j < fake_blocks_out_ids[i].size(); j++) {
            blocks_cells_out_lens[i].push_back(blocks_cells_out_idxs[i][j].size());
        }
    }

    blocks_cells_in_lens = vector<vector<int> >(mpiInfo.comm_size); // длина для каждого процесса -- количсевто блоков от него
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        blocks_cells_in_lens[i] = vector<int>(fake_blocks_in_ids[i].size());
    }

    MPI_Request send_reqs[active_neighs_num];
    MPI_Status send_statuses[active_neighs_num];
    int req_num = 0;
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ( (n != mpiInfo.comm_rank) && (fake_blocks_out_ids[n].size() > 0) ) {
            MPI_Isend(&blocks_cells_out_lens[n][0], fake_blocks_out_ids[n].size(), MPI_INT, n, 0, mpiInfo.comm, &send_reqs[req_num]);
            req_num++;
        }
    }
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ( (n != mpiInfo.comm_rank) && (fake_blocks_in_ids[n].size() > 0) ) {
            MPI_Status status;
            MPI_Recv(&blocks_cells_in_lens[n][0], fake_blocks_in_ids[n].size(), MPI_INT, n, MPI_ANY_TAG, mpiInfo.comm, &status);
        }
    }
    MPI_Waitall(active_neighs_num, send_reqs, send_statuses);

    // формирование исходящих и входящих структур для списков индексов по блокам
    cells_in_idxs = vector<vector<GlobalNumber_t > >(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        int sum = 0;
        for (int j = 0; j < blocks_cells_in_lens[i].size(); j++) {
            sum += blocks_cells_in_lens[i][j];
        }
        cells_in_idxs[i] = vector<GlobalNumber_t>(sum);
    }

    cells_out_idxs = vector<vector<GlobalNumber_t > >(mpiInfo.comm_size); // по сути линеаризованный для каждого процесса blocks_cells_out_idxs
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        for (int j = 0; j < fake_blocks_out_ids[i].size(); j++) {
            GlobalNumber_t blk_j = fake_blocks_out_ids[i][j];
            for (int k = 0; k < blocks_cells_out_lens[i][j]; j++) {
                cells_out_idxs[i].push_back(blocks_cells_out_idxs[i][blk_j][k]);
            }
        }
    }

    // обмен индексами ячеек
    req_num = 0;
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if (fake_blocks_out_ids[n].size() > 0) {
            MPI_Isend(&cells_out_idxs[n][0], cells_out_idxs[n].size(), MPI_LONG_LONG_INT, n, 0, mpiInfo.comm, &send_reqs[req_num]);
            req_num++;
        }
    }
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ( (n != mpiInfo.comm_rank) && (fake_blocks_in_ids[n].size() > 0) ) {
            MPI_Status status;
            MPI_Recv(&cells_in_idxs[n][0], cells_in_idxs[n].size(), MPI_LONG_LONG_INT, n, MPI_ANY_TAG, mpiInfo.comm, &status);
        }
    }
    MPI_Waitall(active_neighs_num, send_reqs, &send_statuses[0]);


    // (3) построение структур для значений температур
    cells_in_idxs_temps = vector<vector<double> >(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        cells_in_idxs_temps[i] = vector<double>(cells_in_idxs[i].size());
    }

    cells_out_idxs_temps = vector<vector<double> >(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        cells_out_idxs_temps[i] = vector<double>(cells_out_idxs[i].size());
    }

}

void Proc::build_needed_cells_for_blocks_map() {
    blocks_cells_out_idxs = vector<map<GlobalNumber_t, vector<GlobalNumber_t > > >(mpiInfo.comm_size);
    for (int blk_i = 0; blk_i < mesh.blocks.size(); blk_i++) {
        int sz = mesh.blocks[blk_i].sz;
        int c_lvl = mesh.blocks[blk_i].cells_lvl;

        // bottom border
        for (int j = 0; j < sz; j++) {
            GlobalNumber_t c_glob_idx = get_cell_global_number_in_full_grid(mesh.blocks[blk_i].idx.get_global_number(), 0, j, c_lvl);

            for (int jjj = 0; jjj < mesh.blocks[blk_i].neighs_down_idxs.size(); jjj++) {
                GlobalNumber_t neigh_blk_i = mesh.blocks[blk_i].neighs_down_idxs[jjj];
                int o = find_owner(neigh_blk_i);
                if (o != mpiInfo.comm_rank) {
                    BlockOfCells* neigh_blk = fake_ghost_blocks[neigh_blk_i];
                    vector<GlobalNumber_t> c_neighs = find_cell_neighs_ids_in_blk(c_glob_idx, c_lvl, neigh_blk, DOWN);
                    
                    for (int jjj2 = 0; jjj2 < c_neighs.size(); jjj2++) {
                        GlobalNumber_t cc = c_neighs[jjj2];
                        if (blocks_cells_out_idxs[o].find(neigh_blk_i) == blocks_cells_out_idxs[o].end()) {
                            blocks_cells_out_idxs[o][neigh_blk_i] = vector<GlobalNumber_t>();
                        }
                        blocks_cells_out_idxs[o][neigh_blk_i].push_back(cc);
                    }
                }
            }
        }
        // upper border
        for (int j = 0; j < sz; j++) {
            GlobalNumber_t c_glob_idx = get_cell_global_number_in_full_grid(mesh.blocks[blk_i].idx.get_global_number(), (sz-1), j, c_lvl);

            for (int jjj = 0; jjj < mesh.blocks[blk_i].neighs_upper_idxs.size(); jjj++) {
                GlobalNumber_t n_blk_i = mesh.blocks[blk_i].neighs_upper_idxs[jjj];
                int o = find_owner(n_blk_i);
                if (o != mpiInfo.comm_rank) {
                    BlockOfCells* n_blk = fake_ghost_blocks[n_blk_i];
                    vector<GlobalNumber_t> c_neighs = find_cell_neighs_ids_in_blk(c_glob_idx, c_lvl, n_blk, UP);
                    for (int jjj2 = 0; jjj2 < c_neighs.size(); jjj2++) {
                        GlobalNumber_t cc = c_neighs[jjj2];
                        if (blocks_cells_out_idxs[o].find(n_blk_i) == blocks_cells_out_idxs[o].end()) {
                            blocks_cells_out_idxs[o][n_blk_i] = vector<GlobalNumber_t>();
                        }
                        blocks_cells_out_idxs[o][n_blk_i].push_back(cc);
                    }
                }
            }
        }
        // left border
        for (int i = 0; i < sz; i++) {
            GlobalNumber_t c_glob_idx = get_cell_global_number_in_full_grid(mesh.blocks[blk_i].idx.get_global_number(), i, 0, c_lvl);

            for (int jjj = 0; jjj < mesh.blocks[blk_i].neighs_left_idxs.size(); jjj++) {
                GlobalNumber_t n_blk_i = mesh.blocks[blk_i].neighs_left_idxs[jjj];
                int o = find_owner(n_blk_i);
                if (o != mpiInfo.comm_rank) {
                    BlockOfCells* n_blk = fake_ghost_blocks[n_blk_i];
                    vector<GlobalNumber_t> c_neighs = find_cell_neighs_ids_in_blk(c_glob_idx, c_lvl, n_blk, LEFT);
                    for (int jjj2 = 0; jjj2 < c_neighs.size(); jjj2++) {
                        GlobalNumber_t cc = c_neighs[jjj2];
                        if (blocks_cells_out_idxs[o].find(n_blk_i) == blocks_cells_out_idxs[o].end()) {
                            blocks_cells_out_idxs[o][n_blk_i] = vector<GlobalNumber_t>();
                        }
                        blocks_cells_out_idxs[o][n_blk_i].push_back(cc);
                    }
                }
            }
        }
        // right border
        for (int i = 0; i < sz; i++) {
            GlobalNumber_t c_glob_idx = get_cell_global_number_in_full_grid(mesh.blocks[blk_i].idx.get_global_number(), i, (sz-1), c_lvl);

            for (int jjj = 0; jjj < mesh.blocks[blk_i].neighs_right_idxs.size(); jjj++) {
                GlobalNumber_t n_blk_i = mesh.blocks[blk_i].neighs_right_idxs[jjj];
                int o = find_owner(n_blk_i);
                if (o != mpiInfo.comm_rank) {
                    BlockOfCells* n_blk = fake_ghost_blocks[n_blk_i];
                    vector<GlobalNumber_t> c_neighs = find_cell_neighs_ids_in_blk(c_glob_idx, c_lvl, n_blk, RIGHT);
                    for (int jjj2 = 0; jjj2 < c_neighs.size(); jjj2++) {
                        GlobalNumber_t cc = c_neighs[jjj2];
                        if (blocks_cells_out_idxs[o].find(n_blk_i) == blocks_cells_out_idxs[o].end()) {
                            blocks_cells_out_idxs[o][n_blk_i] = vector<GlobalNumber_t>();
                        }
                        blocks_cells_out_idxs[o][n_blk_i].push_back(cc);
                    }
                }
            }
        }
    }
}



int Proc::find_owner(GlobalNumber_t cell_id) {
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        MetaInfo proc_meta = meta.GetMetaOfProc(i);
        if ((proc_meta.procStart <= cell_id) && (proc_meta.procEnd >= cell_id)) {
            return i;
        }
    }
    
    std::cout << "owner not found for cell_id=" << cell_id << std::endl;
    return -1;
}


void Proc::build_border_cells_neighs_ptrs() {

}

int Proc::StartExchangeGhosts() {
    
    // cout << mpiInfo.comm_rank << " ExchangeGhosts started\n";
    stat.timers["exchange_ghosts"].Start();

    int temp_l_corr = (time_step_n-1) % 2; // чтобы брать значение temp[0] или temp[1] (-1 т.к. увеличиваем в самом начале шага)
    int cur_temp_idx = temp_l_corr;
    int req_num;

    // non-blocking send to all neighs

    req_num = 0;
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ((n == mpiInfo.comm_rank) || (fake_blocks_out_ids[n].size() == 0)) {
            continue;
        }

        int offset_sum = 0;
        for (int j = 0; j < fake_blocks_out_ids[n].size(); j++) {
            BlockOfCells *blk = mesh.find_block(fake_blocks_out_ids[n][j]);
            for (int k = 0; k < blocks_cells_out_lens[n][j]; k++) {
                GlobalNumber_t cell_idx = cells_out_idxs[n][offset_sum+k];
                SimpleCell *c = blk->find_border_cell_by_global_idx(cell_idx);
                cells_out_idxs_temps[n][offset_sum+k] = c->temp[cur_temp_idx];
            }
            offset_sum += blocks_cells_out_lens[n][j];
        }

        MPI_Isend(&cells_out_idxs_temps[n][0], cells_out_idxs_temps[n].size(), MPI_DOUBLE, n, 0, mpiInfo.comm, &send_reqs[req_num]);
        MPI_Irecv(&cells_in_idxs_temps[n][0], cells_in_idxs_temps[n].size(), MPI_DOUBLE, n, 0, mpiInfo.comm, &recv_reqs[req_num]);

        req_num++;
    }

    // if (FULL_DEBUG) {
    //     cout << mpiInfo.comm_rank << " GHOSTS={ ";
    //     for (int n = 0; n < mpiInfo.comm_size; n++) {
    //         cout << "[ ";
    //         for (int i = 0; i < ghosts_in[n].cells.size(); i++) {
    //             cout << "(" << ghosts_in[n].cells[i].lvl << "," << ghosts_in[n].cells[i].i << "," << ghosts_in[n].cells[i].j << ": " <<  ghosts_in[n].cells[i].temp[cur_temp_idx] << ") ";
    //         }
    //         cout << "] ";
    //     }
    //     cout << "}\n";
    // }

    stat.timers["exchange_ghosts"].Stop();
    // cout << mpiInfo.comm_rank << " ExchangeGhosts finished\n";
    return 0;
}

int Proc::StopExchangeGhosts() {
    stat.timers["exchange_ghosts"].Start();

    int temp_l_corr = (time_step_n-1) % 2; // чтобы брать значение temp[0] или temp[1] (потому что эта функция вызывается уже после обновления номера шага)
    int cur_temp_idx = temp_l_corr;
    int next_temp_idx = (temp_l_corr + 1) % 2;

    MPI_Waitall(active_neighs_num, &send_reqs[0], &send_statuses[0]);
    MPI_Waitall(active_neighs_num, &recv_reqs[0], &recv_statuses[0]);

    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ((n == mpiInfo.comm_rank) || (fake_blocks_in_ids[n].size() == 0)) {
            continue;
        }

        int offset_sum = 0;
        for (int j = 0; j < fake_blocks_in_ids[n].size(); j++) {
            BlockOfCells *blk = mesh.find_block(fake_blocks_in_ids[n][j]);
            for (int k = 0; k < blocks_cells_in_lens[n][j]; k++) {
                GlobalNumber_t cell_idx = cells_in_idxs[n][offset_sum + k];
                SimpleCell *c = (*blk).find_border_cell_by_global_idx(cell_idx);
                c->temp[cur_temp_idx] = cells_in_idxs_temps[n][offset_sum + k];
            }
            offset_sum += blocks_cells_in_lens[n][j];
        }
    }

    stat.timers["exchange_ghosts"].Stop();
    return 0;
}


void Proc::MakeStep() {
    time_step_n++;
    // usleep(10000000);

    stat.timers["step"].Start();
    int temp_l_corr = (time_step_n-1) % 2; // чтобы брать значение temp[0] или temp[1]
    int cur_temp_idx = temp_l_corr;
    int next_temp_idx = (temp_l_corr + 1) % 2;

    StartExchangeGhosts();

    // PrintMyCells();
    // PrintGhostCells();

    // внутренние ячейки блоков
    stat.timers["compute_temps"].Start();
    process_blocks_inner_cells();
    stat.timers["compute_temps"].Stop();

    StopExchangeGhosts();

    // границы блоков
    stat.timers["compute_temps"].Start();
    stat.timers["compute_temps_border"].Start();
    process_blocks_border_cells();
    stat.timers["compute_temps_border"].Stop();
    stat.timers["compute_temps"].Stop();

    stat.timers["step"].Stop();
}

void Proc::process_blocks_inner_cells() {
    int temp_l_corr = (time_step_n-1) % 2; // -1 т.к. увеличиваем счетчик time_step_n в самом начале выполнения шага
    int cur_temp_idx = temp_l_corr;
    int next_temp_idx = (temp_l_corr + 1) % 2;
    
    int blk_i, i, j;

    // внутренние ячейки блоков
    #ifdef USE_OPEN_MP
    # pragma omp parallel
    #endif
    {
    
        int thread_num = 0;
        #ifdef USE_OPEN_MP
        thread_num = omp_get_thread_num();
        #endif
        // cout << "proc " << mpiInfo.comm_rank << " thread " << thread_num << endl;

        #ifdef USE_OPEN_MP
        # pragma omp for private(blk_i, i, j)
        #endif
        for (blk_i = 0; blk_i < mesh.blocks.size(); blk_i++) {

            double d = get_lvl_dx(mesh.blocks[blk_i].cells_lvl);
            int sz = mesh.blocks[blk_i].sz;
            double l = d;
            double S = d * d;

            for (i = 1; i < mesh.blocks[blk_i].sz - 1; i++) {
                for (j = 1; j < mesh.blocks[blk_i].sz - 1; j++) {
                    double flows_sum = 0.0;
                    int n_neighs = 4;
                    double t0 = mesh.blocks[blk_i].cells[i*sz+j].temp[cur_temp_idx];

                    #ifdef USE_SIMPLE_FLOWS
                        flows_sum += mesh.blocks[blk_i].cells[i*sz+(j-1)].temp[cur_temp_idx];
                        flows_sum += mesh.blocks[blk_i].cells[i*sz+(j+1)].temp[cur_temp_idx];
                        flows_sum += mesh.blocks[blk_i].cells[(i-1)*sz+j].temp[cur_temp_idx];
                        flows_sum += mesh.blocks[blk_i].cells[(i+1)*sz+j].temp[cur_temp_idx];
                    #else
                        double t1 = mesh.blocks[blk_i].cells[i*sz+(j-1)].temp[cur_temp_idx];
                        flows_sum += - Area::a * (t1 - t0) / d * l;
                        t1 = mesh.blocks[blk_i].cells[i*sz+(j+1)].temp[cur_temp_idx];
                        flows_sum += - Area::a * (t1 - t0) / d * l;
                        t1 = mesh.blocks[blk_i].cells[(i-1)*sz+j].temp[cur_temp_idx];
                        flows_sum += - Area::a * (t1 - t0) / d * l;
                        t1 = mesh.blocks[blk_i].cells[(i+1)*sz+j].temp[cur_temp_idx];
                        flows_sum += - Area::a * (t1 - t0) / d * l;
                    #endif

                    double x, y;
                    mesh.blocks[blk_i].get_spacial_coords(i, j, &x, &y);
                    double q = Area::Q(x, y, tau * time_step_n);
                    
                    // cout << "flows_sum=" << flows_sum << " q=" << q << "t_new=" << t0 + tau * (flows_sum + q) / S;
                    // cout << "internal block cells: t0=" << t0 << " flows_sum=" << flows_sum << " n_neighs=" << n_neighs << " q=" << q << endl;
                    #ifdef USE_SIMPLE_FLOWS
                        // просто среднее значение соседних ячеек + предыдущее значение в этой + источник
                        flows_sum *= Area::a;
                        mesh.blocks[blk_i].cells[i*sz+j].temp[next_temp_idx] = t0 + flows_sum / n_neighs + q;
                    #else
                        mesh.blocks[blk_i].cells[i*sz+j].temp[next_temp_idx] = t0 + tau * (flows_sum + q) / S;
                    #endif
                }
            }
        }
    }
}

void Proc::process_blocks_border_cells() {

    int temp_l_corr = (time_step_n-1) % 2; // чтобы брать значение temp[0] или temp[1]
    int cur_temp_idx = temp_l_corr;
    int next_temp_idx = (temp_l_corr + 1) % 2;

    int i, j;
    for (int blk_i = 0; blk_i < mesh.blocks.size(); blk_i++) {
        BlockOfCells& blk = mesh.blocks[blk_i];

        // upper border
        i = blk.sz - 1;
        for (j = 0; j < blk.sz ; j++) {
            process_block_border_cell(blk, i, j);
        }
        // bottom border
        i = 0;
        for (j = 0; j < blk.sz ; j++) {
            process_block_border_cell(blk, i, j);
        }
        // right border
        j = 0;
        for (i = 1; i < blk.sz - 1 ; i++) {
            process_block_border_cell(blk, i, j);
        }
         // left border
        j = blk.sz - 1;
        for (i = 1; i < blk.sz - 1 ; i++) {
            process_block_border_cell(blk, i, j);
        }
    }
}

 bool is_cell_on_area_border(BlockOfCells &blk, int i, int j, Neigh dir) {
    switch (dir) {
    case DOWN:
        return (i == 0) && (blk.neighs_down_idxs.size() == 0);
    case UP:
        return (i == blk.sz-1) && (blk.neighs_upper_idxs.size() == 0);
    case LEFT:
        return (j == 0) && (blk.neighs_left_idxs.size() == 0);
    case RIGHT:
        return (j == blk.sz-1) && (blk.neighs_right_idxs.size() == 0);
    default:
        return false;
    }
}

std::pair<double,int> Proc::compute_inblock_flows_for_border_cell(BlockOfCells& blk, int i, int j) {
    int cur_temp_idx = (time_step_n-1) % 2;
    int sz = blk.sz;

    double d = get_lvl_dx(blk.cells_lvl);
    double l = d;

    double t0 = blk.cells[i*sz+j].temp[cur_temp_idx];

    double flows_this_block_sum = 0.0;
    int n_neighs = 0;

    if (i > 0) { // down
        n_neighs++;
        double t1 = blk.cells[(i-1)*sz+j].temp[cur_temp_idx];
        #ifdef USE_SIMPLE_FLOWS
            flows_this_block_sum += t1;
        #else
            flows_this_block_sum += - Area::a * (t1 - t0) / d * l;
        #endif
    }
    if (i < sz-1) { // up
        n_neighs++;
        double t1 = blk.cells[(i+1)*sz+j].temp[cur_temp_idx];
        #ifdef USE_SIMPLE_FLOWS
            flows_this_block_sum += t1;
        #else
            flows_this_block_sum += - Area::a * (t1 - t0) / d * l;
        #endif
    }
    if (j > 0) { // left
        n_neighs++;
        double t1 = blk.cells[i*sz+(j-1)].temp[cur_temp_idx];
        #ifdef USE_SIMPLE_FLOWS
            flows_this_block_sum += t1;
        #else
            flows_this_block_sum += - Area::a * (t1 - t0) / d * l;
        #endif
    }
    if (j < sz - 1) { // right
        n_neighs++;
        double t1 = blk.cells[i*sz+(j+1)].temp[cur_temp_idx];
        #ifdef USE_SIMPLE_FLOWS
            flows_this_block_sum += t1;
        #else
            flows_this_block_sum += - Area::a * (t1 - t0) / d * l;
        #endif
    }

    return std::pair<double,int>(flows_this_block_sum, n_neighs);
}

std::pair<double,int> Proc::compute_border_flows_for_border_cell(BlockOfCells& blk, int i, int j) {
    int cur_temp_idx = (time_step_n-1) % 2;
    double t0 = blk.cells[i*blk.sz+j].temp[cur_temp_idx];

    double x, y;
    blk.get_spacial_coords(i, j, &x, &y);

    double d = get_lvl_dx(blk.cells_lvl);
    double l = d;
    
    double flows_borders_sum = 0.0;
    int n_neighs = 0;

    char border_cond_type;
    double (*cond_func)(double, double, double);
    for (int border_dir = DOWN; border_dir != END; border_dir++) {
        if (is_cell_on_area_border(blk, i, j, Neigh(border_dir))) {
            stat.timers["get_border_cond"].Start();
            get_border_cond(&border_cond_type, &cond_func, Neigh(border_dir));
            stat.timers["get_border_cond"].Stop();

            double t1 = cond_func(x+d/2, y, time_step_n * tau);
            if (border_cond_type == 1) {
                #ifdef USE_SIMPLE_FLOWS
                    flows_borders_sum += t1;
                    n_neighs++;
                #else
                    flows_borders_sum += - Area::a * (t1 - t0) / (d/2) * l;
                #endif
            } 
            // else if (border_cond_type == 2) {
            //     flows_borders_sum += t1;
            // }
        }
    }

    return std::pair<double,int>(flows_borders_sum, n_neighs);
}

vector<Neigh> Proc::get_possible_adjacent_blocks_dirs_for_cell(BlockOfCells& blk, int i, int j) {
    vector<Neigh> dirs;
    if (i == 0) {
        dirs.push_back(DOWN);
    }
    if (i == blk.sz - 1) {
        dirs.push_back(UP);
    }
    if (j == 0) {
        dirs.push_back(LEFT);
    }
    if (j == blk.sz - 1) {
        dirs.push_back(RIGHT);
    }

    return dirs;
}

vector<GlobalNumber_t> get_possible_adjacent_blocks(BlockOfCells& blk, Neigh neigh_dir) {
    switch(neigh_dir) {
    case LEFT:
        return blk.neighs_left_idxs;
    case RIGHT:
        return blk.neighs_right_idxs;
    case UP:
        return blk.neighs_upper_idxs;
    case DOWN:
        return blk.neighs_down_idxs;
    default:
        return vector<GlobalNumber_t>();
    }
}

std::pair<double,int> Proc::compute_interblock_flows_for_border_cell(BlockOfCells& blk, int i, int j) {
    double flows_other_blocks_sum = 0.0;
    int n_neighs = 0;

    vector<Neigh> dirs = get_possible_adjacent_blocks_dirs_for_cell(blk, i, j);
    for (std::vector<Neigh>::iterator it = dirs.begin(); it != dirs.end(); it++) {
        // if (dir == DOWN) {
        //     cout << "computing DOWN border interblock flows ("<< i << "," << j<<") neighs_down_sz=" << blk.neighs_down_idxs.size() << endl;
        // }
        std::pair<double,int> tmp_pair = compute_interblock_flows_for_border_cell(blk, i, j, *it);
        flows_other_blocks_sum += tmp_pair.first;
        n_neighs += tmp_pair.second;
    }

    return std::pair<double,int>(flows_other_blocks_sum, n_neighs);
} 

std::pair<double,int> Proc::compute_interblock_flows_for_border_cell(BlockOfCells& blk, int i, int j, Neigh neigh_dir) {
    int cur_temp_idx = (time_step_n-1) % 2;
    int sz = blk.sz;

    double d = get_lvl_dx(blk.cells_lvl);
    double l = d;

    double t0 = blk.cells[i*sz+j].temp[cur_temp_idx];


    double flows_other_blocks_sum = 0.0;
    int n_neighs = 0;

    GlobalNumber_t c_glob_idx = get_cell_global_number_in_full_grid(blk.idx.get_global_number(), i, j, blk.cells_lvl);
    vector<GlobalNumber_t> n_blks_is = get_possible_adjacent_blocks(blk, neigh_dir);
    for (int nb_i = 0; nb_i < n_blks_is.size(); nb_i++) {
        GlobalNumber_t n_blk_i = n_blks_is[nb_i];
        int n_blk_owner = find_owner(n_blk_i);
        bool is_my_blk = (n_blk_owner == mpiInfo.comm_rank);
        BlockOfCells *n_blk = is_my_blk ? mesh.find_block(n_blk_i) : fake_ghost_blocks[n_blk_i];
        vector<GlobalNumber_t> cell_neighs = find_cell_neighs_ids_in_blk(c_glob_idx, blk.cells_lvl, n_blk, neigh_dir);
        // if (neigh_dir == DOWN) {
        //     cout << "DOWN cell_neighs: ";
        //     for (auto gn : cell_neighs) {
        //         cout << gn << " ";
        //     }
        //     cout << endl;
        // }
        for (int nc_i = 0; nc_i < cell_neighs.size(); nc_i++) {
            GlobalNumber_t cc = cell_neighs[nc_i];
            double t1 = 0.0;
            // если уровень ячеек в этом соседнем блоке больше, то ячейки там меньше, то длина грани -- половинка
            if (n_blk->cells_lvl > blk.cells_lvl) {
                l = l / 2;
            }

            bool found = false;
            if (is_my_blk) {
                SimpleCell *c = n_blk->find_border_cell_by_global_idx(cc);
                t1 = c->temp[cur_temp_idx];
                found = true;
            } else {
                for (int k = 0; k < cells_in_idxs[n_blk_owner].size(); k++) {
                    if (cells_in_idxs[n_blk_owner][k] == cc) {
                        t1 = cells_in_idxs_temps[n_blk_owner][k];
                        found = true;
                    }
                }
            }
            if (found) {
                #ifdef USE_SIMPLE_FLOWS
                    flows_other_blocks_sum += t1;
                    n_neighs++;
                #else
                    flows_other_blocks_sum += - Area::a * (t1 - t0) / d * l;
                #endif
            }
        }
    }

    // if (neigh_dir == DOWN) {
    //     cout << "DOWNflow=" << flows_other_blocks_sum << endl;
    // }
    // if (neigh_dir == LEFT) {
    //     cout << "LEFTflow=" << flows_other_blocks_sum << endl;
    // }

    return std::pair<double,int>(flows_other_blocks_sum, n_neighs);
}

void Proc::process_block_border_cell(BlockOfCells& blk, int i, int j) {
    int cur_temp_idx = (time_step_n-1) % 2;
    int next_temp_idx = (cur_temp_idx + 1) % 2;

    double d = get_lvl_dx(blk.cells_lvl);
    double S = d * d;

    double flows_sum = 0.0;
    int n_neighs = 0;

    // Потоки: 
    //   1) от соседних ячеек в этом же блоке
    //   2) от границ области
    //   3) от соседних ячеек в других блоках (на этом же проце или нет)

    // neighs in this block
    // stat.timers["compute_temps_border_inblock"].Start();
    std::pair<double, int> tmp_pair = compute_inblock_flows_for_border_cell(blk, i, j);
    flows_sum += tmp_pair.first;
    n_neighs += tmp_pair.second;
    // stat.timers["compute_temps_border_inblock"].Stop();
    

    // border flows
    // stat.timers["compute_temps_border_edges"].Start();
    tmp_pair = compute_border_flows_for_border_cell(blk, i, j);
    flows_sum += tmp_pair.first;
    n_neighs += tmp_pair.second;
    // stat.timers["compute_temps_border_edges"].Stop();

    // neighs in other blocks
    // stat.timers["compute_temps_border_interblock"].Start();
    tmp_pair = compute_interblock_flows_for_border_cell(blk, i, j);
    flows_sum += tmp_pair.first;
    n_neighs += tmp_pair.second;
    // stat.timers["compute_temps_border_interblock"].Stop();


    // result
    // stat.timers["compute_temps_border_res"].Start();
    double x, y;
    blk.get_spacial_coords(i, j, &x, &y);
    double q = Area::Q(x, y, tau * time_step_n);

    int c_blk_num = i*blk.sz+j;
    double t0 = blk.cells[c_blk_num].temp[cur_temp_idx];
    #ifdef USE_SIMPLE_FLOWS
        flows_sum *= Area::a;
        blk.cells[c_blk_num].temp[next_temp_idx] = t0 + flows_sum / n_neighs + q;
    #else
        blk.cells[c_blk_num].temp[next_temp_idx] = t0 + tau * (flows_sum + q) / S;
    #endif
    // stat.timers["compute_temps_border_res"].Stop();
}



void Proc::get_border_cond(char *cond_type, double (**cond_func)(double, double, double), Neigh border) {
    if (border == DOWN) {
        Area::get_border_cond(Area::DOWN, cond_type, cond_func);
    } else if (border == UP) {
        Area::get_border_cond(Area::UP, cond_type, cond_func);
    } else if (border == RIGHT) {
        Area::get_border_cond(Area::RIGHT, cond_type, cond_func);
    } else if (border == LEFT) {
        Area::get_border_cond(Area::LEFT, cond_type, cond_func);
    }
}

void Proc::WriteT(char * filename) {
    stat.timers["io"].Start();

    vector<char> buf = mesh.GenWriteStruct(1);
    int len = buf.size();
    int *lens = new int[mpiInfo.comm_size];
//    stat.timers["communication"].Start();
    MPI_Allgather(&len, 1, MPI_INT, (void *)lens, 1, MPI_INT, mpiInfo.comm);
//    stat.timers["communication"].Stop();
    int offset = 0;
    for (int i = 0; i < mpiInfo.comm_rank; i++) {
        offset += lens[i];
    }

    // std::cout << mpiInfo.comm_rank << " offset=" << offset << " wr_len=" << len << std::endl;

    MPI_File fh;
    MPI_File_open( mpiInfo.comm, 
                filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_at(fh, offset, &buf[0], len, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    stat.timers["io"].Stop();
}

void Proc::WriteStat() {
    std::cout << "MpiInfo: " << mpiInfo.toString() << std::endl;
    double step = stat.timers["step"].FullDur();
    double compute_temps = stat.timers["compute_temps"].FullDur();
    double compute_temps_border = stat.timers["compute_temps_border"].FullDur();
    double compute_temps_border_inblock = stat.timers["compute_temps_border_inblock"].FullDur();
    double compute_temps_border_edges = stat.timers["compute_temps_border_edges"].FullDur();
    double compute_temps_border_interblock = stat.timers["compute_temps_border_interblock"].FullDur();
    double compute_temps_border_res = stat.timers["compute_temps_border_res"].FullDur();
    double exchange_ghosts = stat.timers["exchange_ghosts"].FullDur();
    double build_ghosts = stat.timers["build_ghosts"].FullDur();
    double init_mesh = stat.timers["init_mesh"].FullDur();

    double max_step, 
        max_compute_temps, 
        max_compute_temps_border, 
        max_compute_temps_border_inblock, 
        max_compute_temps_border_edges, 
        max_compute_temps_border_interblock, 
        max_compute_temps_border_res, 
        max_exchange_ghosts, 
        max_build_ghosts, 
        max_init_mesh;

    MPI_Reduce(&step, &max_step, 1, MPI_DOUBLE, MPI_MAX, 0, mpiInfo.comm);
    MPI_Reduce(&compute_temps, &max_compute_temps, 1, MPI_DOUBLE, MPI_MAX, 0, mpiInfo.comm);
    MPI_Reduce(&compute_temps_border, &max_compute_temps_border, 1, MPI_DOUBLE, MPI_MAX, 0, mpiInfo.comm);
    MPI_Reduce(&compute_temps_border_inblock, &max_compute_temps_border_inblock, 1, MPI_DOUBLE, MPI_MAX, 0, mpiInfo.comm);
    MPI_Reduce(&compute_temps_border_edges, &max_compute_temps_border_edges, 1, MPI_DOUBLE, MPI_MAX, 0, mpiInfo.comm);
    MPI_Reduce(&compute_temps_border_interblock, &max_compute_temps_border_interblock, 1, MPI_DOUBLE, MPI_MAX, 0, mpiInfo.comm);
    MPI_Reduce(&compute_temps_border_res, &max_compute_temps_border_res, 1, MPI_DOUBLE, MPI_MAX, 0, mpiInfo.comm);
    MPI_Reduce(&exchange_ghosts, &max_exchange_ghosts, 1, MPI_DOUBLE, MPI_MAX, 0, mpiInfo.comm);
    MPI_Reduce(&build_ghosts, &max_build_ghosts, 1, MPI_DOUBLE, MPI_MAX, 0, mpiInfo.comm);
    MPI_Reduce(&init_mesh, &max_init_mesh, 1, MPI_DOUBLE, MPI_MAX, 0, mpiInfo.comm);

    if (mpiInfo.comm_rank == 0) {
        cout << "\n\n\n\nSTAT| max_step: " << max_step
                << " | max_compute_temps: " << max_compute_temps
                << " | max_compute_temps_border: " << max_compute_temps_border
                << " | max_compute_temps_border_inblock: " << max_compute_temps_border_inblock
                << " | max_compute_temps_border_edges: " << max_compute_temps_border_edges
                << " | max_compute_temps_border_interblock: " << max_compute_temps_border_interblock
                << " | max_compute_temps_border_res: " << max_compute_temps_border_res
                << " | max_exchange_ghosts: " << max_exchange_ghosts
                << " | max_build_ghosts: " << max_build_ghosts
                << " | max_init_mesh: " << max_init_mesh << " |\n\n\n\n" << endl; 
    }

}


void Proc::PrintMyBlocks() {

    int temp_l_corr = time_step_n % 2; // чтобы брать значение temp[0] или temp[1]
    int cur_temp_idx = temp_l_corr;

    // cout << mpiInfo.comm_rank << " CELLS={ ";
    // for (int i = 0; i < mesh.blocks.size(); i++) {
    //     cout << "(" << mesh.blocks[i].idx.lvl << "," << mesh.blocks[i].idx.i << "," << mesh.blocks[i].idx.j <<
    //           mesh.blocks[i].cells_lvl << ", " << mesh.blocks[i].sz << ") ";
    // }
    // cout << "}\n";
}

void Proc::PrintGhostCells() {

//    int temp_l_corr = time_step_n % 2; // чтобы брать значение temp[0] или temp[1]
//    int cur_temp_idx = temp_l_corr;
//
//    cout << mpiInfo.comm_rank << " GHOSTS={ ";
//        for (int n = 0; n < mpiInfo.comm_size; n++) {
//            cout << n << ":[ ";
//            for (int i = 0; i < ghosts_in[n].cells.size(); i++) {
//                cout << "(" << ghosts_in[n].cells[i].lvl << "," << ghosts_in[n].cells[i].i << "," << ghosts_in[n].cells[i].j << ": " <<  ghosts_in[n].cells[i].temp[cur_temp_idx] << ") ";
//            }
//            cout << "] ";
//        }
//        cout << "}\n";
}

size_t Proc::GetProcAllocMem() {
    size_t res = 0;
    for (int i = 0; i < mesh.blocks.size(); i++) {
        res += mesh.blocks[i].GetMemSize();
    }

    for (int i = 0; i < mpiInfo.comm_size; i++) {
        res += sizeof(GlobalNumber_t) * fake_blocks_out_ids[i].size();
        res += sizeof(GlobalNumber_t) * fake_blocks_in_ids[i].size();
        res += sizeof(int) * blocks_cells_out_lens[i].size();
        res += sizeof(int) * blocks_cells_in_lens[i].size();
    }

    for (int i = 0; i < mpiInfo.comm_size; i++) {
        res += sizeof(GlobalNumber_t) * cells_out_idxs[i].size();
        res += sizeof(GlobalNumber_t) * cells_in_idxs[i].size();
        res += sizeof(double) * cells_out_idxs_temps[i].size();
        res += sizeof(double) * cells_in_idxs_temps[i].size();
    }

    return res;
}
