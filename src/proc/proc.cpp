#include "proc.h"

#define FULL_DEBUG 0

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
    tau = (min_dx * min_dx) / 4;
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
    for (auto t: timers) {
        stringStream << t.first << ": " << t.second.FullDur() << std::endl;
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
    stat.timers["total"].Start();
    stat.timers["io"] = MpiTimer();
    stat.timers["build_ghosts"] = MpiTimer();
    stat.timers["exchange_ghosts"] = MpiTimer();
    stat.timers["communication"] = MpiTimer();
    stat.timers["find_cell"] = MpiTimer();
    stat.timers["step"] = MpiTimer();
    stat.timers["get_border_cond"] = MpiTimer();
    stat.timers["get_possible_neighs"] = MpiTimer();
    stat.timers["sort_neighs"] = MpiTimer();
}

Proc::~Proc() {

    cout << mpiInfo.comm_rank << " DESTROYING PROC\n";

    stat.timers["total"].Stop();
    std::cout << mpiInfo.comm_rank << " ";
    std::cout << stat.toString() << std::endl;
}

int Proc::MPIInit(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    mpiInfo.comm = MPI_COMM_WORLD;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiInfo.comm_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiInfo.comm_size);

    return 0;
}
int Proc::MPIFinalize() {
    MPI_Finalize();
    return 0;
}

int Proc::InitMesh(string offsets_filename, string blocks_filename, double (*start_func)(double, double)) {
    stat.timers["init_mesh"] = MpiTimer();
    stat.timers["init_mesh"].Start();

    // [0] - offset, [1] - len
    int range[2];
    int one_sz = 5 * sizeof(int);

    MPI_File fh;
    MPI_File_open( mpiInfo.comm, offsets_filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_at(fh, 2 * mpiInfo.comm_rank * sizeof(int), &range, 2, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    std::cout << mpiInfo.comm_rank << " read range[0]=" << range[0] << " range[1]=" << range[1] << std::endl;

    vector<char> buffer(range[1], 1);
    std::cout << "will read cells file\n";

    MPI_File_open( mpiInfo.comm, blocks_filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_read_at(fh, range[0], &buffer[0], range[1]-1, MPI_CHAR, MPI_STATUS_IGNORE);
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

    stat.timers["build_ghosts"].Stop();
    MPI_Barrier(mpiInfo.comm);

    cout << mpiInfo.comm_rank << " built ghosts\n";
    return 0;
}


void Proc::build_fake_ghost_blocks() {

    // строим списки индексов соседей блоков по процессам
    fake_blocks_out_ids = vector<vector<GlobalNumber_t>>(mpiInfo.comm_size);
    for (int blk_i = 0; blk_i < mesh.blocks.size(); blk_i++) {
        for (GlobalNumber_t neight_num: mesh.blocks[blk_i].neighs_left_idxs) {
            int o = find_owner(neight_num);
            if (o != mpiInfo.comm_rank) {
                fake_blocks_out_ids[o].push_back(neight_num);
            }
        }
    }

    // сортируем и удаляем повторения
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        std::sort(fake_blocks_out_ids[i].begin(), fake_blocks_out_ids[i].end());
        auto last = std::unique(fake_blocks_out_ids[i].begin(), fake_blocks_out_ids[i].end());
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
    cout << mpiInfo.comm_rank << " active_neighs_num=" << active_neighs_num << endl;

    send_reqs = vector<MPI_Request>(active_neighs_num);
    send_statuses = vector<MPI_Status>(active_neighs_num);
    recv_reqs = vector<MPI_Request>(active_neighs_num);
    recv_statuses = vector<MPI_Status>(active_neighs_num);


    // обмениваемся длинами чтобы знать сколько блоков от кого принимать
    vector<int> in_lens(mpiInfo.comm_size, 0);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        MPI_Gather(&out_lens[i], 1, MPI_INT, &in_lens[0], 1, MPI_INT, i, mpiInfo.comm);
    }


    // (1) обмен айдишниками блоков
    // формируем массивы для приема
    fake_blocks_in_ids = vector<vector<GlobalNumber_t>>(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        fake_blocks_in_ids[i] = vector<GlobalNumber_t>(in_lens[i]);
    }

    // отправляем и принимаем
    int req_num = 0;
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if (out_lens[n] > 0) {
            MPI_Isend(&fake_blocks_out_ids[n][0], out_lens[n], MPI_LONG_LONG_INT, n, 0, mpiInfo.comm, &send_reqs[req_num]);
            req_num++;
        }
    }
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ( (n != mpiInfo.comm_rank) && (in_lens[n] > 0) ) {
            MPI_Status status;
            MPI_Recv(&fake_blocks_in_ids[n][0], in_lens[n], MPI_LONG_LONG_INT, n, MPI_ANY_TAG, mpiInfo.comm, &status);
        }
    }
    MPI_Waitall(active_neighs_num, &send_reqs[0], &send_statuses[0]);


    // (2) обмен информацией о блоках
    // формируем массивы для отправки
    vector<vector<int>> data_fake_blocks_out;
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        data_fake_blocks_out.push_back(vector<int>());
        for (int j = 0; j < fake_blocks_out_ids[i].size(); j++) {
            BlockOfCells *blk = mesh.find_block(fake_blocks_out_ids[i][j]);
            data_fake_blocks_out[i].push_back(blk->idx.lvl);
            data_fake_blocks_out[i].push_back(blk->cells_lvl);
            data_fake_blocks_out[i].push_back(blk->sz);
        }
    }
    // формируем массивы для приема
    vector<vector<int>> data_fake_blocks_in;
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        data_fake_blocks_in.push_back(vector<int>(in_lens[i]*3));
    }

    // отправляем и принимаем
    req_num = 0;
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if (out_lens[n] > 0) {
            MPI_Isend(&data_fake_blocks_out[n][0], out_lens[n] * 3, MPI_INT, n, 0, mpiInfo.comm, &send_reqs[req_num]);
            req_num++;
        }
    }
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ( (n != mpiInfo.comm_rank) && (in_lens[n] > 0) ) {
            MPI_Status status;
            MPI_Recv(&data_fake_blocks_in[n][0], in_lens[n] * 3, MPI_INT, n, MPI_ANY_TAG, mpiInfo.comm, &status);
        }
    }
    MPI_Waitall(active_neighs_num, &send_reqs[0], &send_statuses[0]);


    // (3) формируем собственно мапу
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        for (int j = 0; j < in_lens[i]; j++) {
            fake_ghost_blocks[fake_blocks_in_ids[i][j]] = new BlockOfCells(
                    data_fake_blocks_in[i][j+1],
                    data_fake_blocks_in[i][j],
                    data_fake_blocks_in[i][j+2],
                    fake_blocks_in_ids[i][j]);
        }
    }
}

void Proc::build_ghost_cells() {

    // (1) построение мапки нужных ячеек по блокам
    build_needed_cells_for_blocks_map();


    // (2) обмен этими нужными ячейками

    // обмен длинами
    blocks_cells_out_lens = vector<vector<int>>(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        for (int j = 0; j < fake_blocks_out_ids[i].size(); j++) {
            blocks_cells_out_lens[i].push_back(cells_out_idxs[i][fake_blocks_out_ids[i][j]].size());
        }
    }

    blocks_cells_in_lens = vector<vector<int>>(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        blocks_cells_in_lens[i] = vector<int>(fake_blocks_in_ids[i].size());
    }

    MPI_Request send_reqs[active_neighs_num];
    MPI_Status send_statuses[active_neighs_num];
    int req_num = 0;
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if (fake_blocks_out_ids[n].size() > 0) {
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
    cells_in_idxs_tmp = vector<vector<GlobalNumber_t >>(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        int sum = 0;
        for (int j = 0; j < blocks_cells_in_lens[i].size(); j++) {
            sum += blocks_cells_in_lens[i][j];
        }
        cells_in_idxs_tmp[i] = vector<GlobalNumber_t>(sum);
    }

    cells_out_idxs_tmp = vector<vector<GlobalNumber_t >>(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        for (int j = 0; j < fake_blocks_out_ids[i].size(); j++) {
            for (int k = 0; k < blocks_cells_out_lens[i][j]; j++) {
                cells_out_idxs_tmp[i].push_back(cells_out_idxs[i][j][k]);
            }
        }
    }

    // обмен индексами ячеек
    req_num = 0;
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if (fake_blocks_out_ids[n].size() > 0) {
            MPI_Isend(&cells_out_idxs_tmp[n][0], cells_out_idxs_tmp[n].size(), MPI_LONG_LONG_INT, n, 0, mpiInfo.comm, &send_reqs[req_num]);
            req_num++;
        }
    }
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ( (n != mpiInfo.comm_rank) && (fake_blocks_in_ids[n].size() > 0) ) {
            MPI_Status status;
            MPI_Recv(&cells_in_idxs_tmp[n][0], cells_in_idxs_tmp[n].size(), MPI_LONG_LONG_INT, n, MPI_ANY_TAG, mpiInfo.comm, &status);
        }
    }
    MPI_Waitall(active_neighs_num, send_reqs, &send_statuses[0]);


    // (3) построение структур для значений температур
    cells_in_idxs_temps = vector<vector<double>>(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        cells_in_idxs_temps[i] = vector<double>(cells_in_idxs_tmp[i].size());
    }

    cells_out_idxs_temps = vector<vector<double>>(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        cells_out_idxs_temps[i] = vector<double>(cells_out_idxs_tmp[i].size());
    }

}

void Proc::build_needed_cells_for_blocks_map() {
    cells_out_idxs = vector<map<GlobalNumber_t, vector<GlobalNumber_t >>>(mpiInfo.comm_size);
    for (int blk_i = 0; blk_i < mesh.blocks.size(); blk_i++) {
        int sz = mesh.blocks[blk_i].sz;
        int c_lvl = mesh.blocks[blk_i].cells_lvl;

        // bottom border
        for (int j = 0; j < sz; j++) {
            GlobalNumber_t c_glob_idx = get_glob_idx(mesh.blocks[blk_i].idx.get_global_number(), 0*sz + j, c_lvl);

            for (GlobalNumber_t n_blk_i: mesh.blocks[blk_i].neighs_down_idxs) {
                int o = find_owner(n_blk_i);
                if (o != mpiInfo.comm_rank) {
                    BlockOfCells* n_blk = fake_ghost_blocks[n_blk_i];
                    vector<GlobalNumber_t> c_neighs = find_cell_neighs_ids_in_blk(c_glob_idx, c_lvl, n_blk, Neigh::DOWN);
                    for (GlobalNumber_t cc: c_neighs) {
                        if (cells_out_idxs[o].find(n_blk_i) == cells_out_idxs[o].end()) {
                            cells_out_idxs[o][n_blk_i] = vector<GlobalNumber_t>();
                        }
                        cells_out_idxs[o][n_blk_i].push_back(cc);
                    }
                }
            }
        }
        // upper border
        for (int j = 0; j < sz; j++) {
            GlobalNumber_t c_glob_idx = get_glob_idx(mesh.blocks[blk_i].idx.get_global_number(), (sz-1)*sz + j, c_lvl);

            for (GlobalNumber_t n_blk_i: mesh.blocks[blk_i].neighs_down_idxs) {
                int o = find_owner(n_blk_i);
                if (o != mpiInfo.comm_rank) {
                    BlockOfCells* n_blk = fake_ghost_blocks[n_blk_i];
                    vector<GlobalNumber_t> c_neighs = find_cell_neighs_ids_in_blk(c_glob_idx, c_lvl, n_blk, Neigh::UP);
                    for (GlobalNumber_t cc: c_neighs) {
                        if (cells_out_idxs[o].find(n_blk_i) == cells_out_idxs[o].end()) {
                            cells_out_idxs[o][n_blk_i] = vector<GlobalNumber_t>();
                        }
                        cells_out_idxs[o][n_blk_i].push_back(cc);
                    }
                }
            }
        }
        // left border
        for (int i = 0; i < sz; i++) {
            GlobalNumber_t c_glob_idx = get_glob_idx(mesh.blocks[blk_i].idx.get_global_number(), i*sz + 0, c_lvl);

            for (GlobalNumber_t n_blk_i: mesh.blocks[blk_i].neighs_down_idxs) {
                int o = find_owner(n_blk_i);
                if (o != mpiInfo.comm_rank) {
                    BlockOfCells* n_blk = fake_ghost_blocks[n_blk_i];
                    vector<GlobalNumber_t> c_neighs = find_cell_neighs_ids_in_blk(c_glob_idx, c_lvl, n_blk, Neigh::DOWN);
                    for (GlobalNumber_t cc: c_neighs) {
                        if (cells_out_idxs[o].find(n_blk_i) == cells_out_idxs[o].end()) {
                            cells_out_idxs[o][n_blk_i] = vector<GlobalNumber_t>();
                        }
                        cells_out_idxs[o][n_blk_i].push_back(cc);
                    }
                }
            }
        }
        // right border
        for (int i = 0; i < sz; i++) {
            GlobalNumber_t c_glob_idx = get_glob_idx(mesh.blocks[blk_i].idx.get_global_number(), i*sz + (sz-1), c_lvl);

            for (GlobalNumber_t n_blk_i: mesh.blocks[blk_i].neighs_down_idxs) {
                int o = find_owner(n_blk_i);
                if (o != mpiInfo.comm_rank) {
                    BlockOfCells* n_blk = fake_ghost_blocks[n_blk_i];
                    vector<GlobalNumber_t> c_neighs = find_cell_neighs_ids_in_blk(c_glob_idx, c_lvl, n_blk, Neigh::DOWN);
                    for (GlobalNumber_t cc: c_neighs) {
                        if (cells_out_idxs[o].find(n_blk_i) == cells_out_idxs[o].end()) {
                            cells_out_idxs[o][n_blk_i] = vector<GlobalNumber_t>();
                        }
                        cells_out_idxs[o][n_blk_i].push_back(cc);
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
    
    // std::cout << "owner not found for cell_id=" << cell_id << std::endl;
    return -1;
}


int Proc::StartExchangeGhosts() {
    
    // cout << mpiInfo.comm_rank << " ExchangeGhosts started\n";
    stat.timers["exchange_ghosts"].Start();

    int temp_l_corr = time_step_n % 2; // чтобы брать значение temp[0] или temp[1]
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
                GlobalNumber_t cell_idx = cells_out_idxs_tmp[n][offset_sum+j];
                SimpleCell *c = (*blk).find_border_cell_by_global_idx(cell_idx);
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

    MPI_Waitall(active_neighs_num, &send_reqs[0], &send_statuses[0]);
    MPI_Waitall(active_neighs_num, &recv_reqs[0], &recv_statuses[0]);

    // todo check statuses



    stat.timers["exchange_ghosts"].Stop();
    return 0;
}


void Proc::MakeStep() {

    // usleep(10000000);

    stat.timers["step"].Start();
    int temp_l_corr = time_step_n % 2; // чтобы брать значение temp[0] или temp[1]
    int cur_temp_idx = temp_l_corr;
    int next_temp_idx = (temp_l_corr + 1) % 2;

    StartExchangeGhosts();

    // PrintMyCells();
    // PrintGhostCells();

    time_step_n++;

    for (int blk_i = 0; blk_i < mesh.blocks.size(); blk_i++) {

        double d = get_lvl_dx(mesh.blocks[blk_i].cells_lvl);
        int sz = mesh.blocks[blk_i].sz;

        // внутренние ячейки блока
        for (int i = 1; i < mesh.blocks[blk_i].sz - 1; i++) {
            for (int j = 1; j < mesh.blocks[blk_i].sz - 1; j++) {
                double flows_sum;
                double t0 = mesh.blocks[blk_i].cells[i*sz+j].temp[cur_temp_idx];

                double t1 = mesh.blocks[blk_i].cells[i*sz+(j-1)].temp[cur_temp_idx];
                flows_sum += - Area::a * (t1 - t0) / d;
                t1 = mesh.blocks[blk_i].cells[i*sz+(j+1)].temp[cur_temp_idx];
                flows_sum += - Area::a * (t1 - t0) / d;
                t1 = mesh.blocks[blk_i].cells[(i-1)*sz+j].temp[cur_temp_idx];
                flows_sum += - Area::a * (t1 - t0) / d;
                t1 = mesh.blocks[blk_i].cells[(i+1)*sz+j].temp[cur_temp_idx];
                flows_sum += - Area::a * (t1 - t0) / d;

                double x, y;
                mesh.blocks[blk_i].get_spacial_coords(i, j, &x, &y);
                double q = Area::Q(x, y, tau * time_step_n);
                double S = (2*d) * (2*d);

                mesh.blocks[blk_i].cells[i*sz+j].temp[next_temp_idx] = t0 - tau * (flows_sum - q) / S;
            }
        }
    }

    StopExchangeGhosts();

    // границы блоков
//    for (int blk_i = 0; blk_i < mesh.blocks.size(); blk_i++){
//
//        double x, y;
//        cell.get_spacial_coords(&x, &y);
//
//        // ищем всех возможных соседей
//        stat.timers["get_possible_neighs"].Start();
//        vector<GlobalNumber_t> possible_neigh_ids = cell.get_all_possible_neighbours_ids();
//        stat.timers["get_possible_neighs"].Stop();
//
//
//        // вычисляем потоки междуу ячейкой и реально существующими соседями
//        vector<double> termo_flows;
//        for (GlobalNumber_t id: possible_neigh_ids) {
//            int owner = find_owner(id);
//            if (owner == -1) {
//                continue;
//            }
//            Cell neigh_cell;
//
//            if (owner == mpiInfo.comm_rank) { // ячейка у меня
//                stat.timers["find_cell"].Start();
//                int idx = mesh.FindCell(id, &neigh_cell);
//                stat.timers["find_cell"].Stop();
//                if (idx == -1) {
//                    continue;
//                }
//            } else { // ячейка-призрак
//                stat.timers["find_cell"].Start();
//                int idx = ghosts_in[owner].FindCell(id, &neigh_cell);
//                stat.timers["find_cell"].Stop();
//                if (idx == -1) {
//                    continue;
//                }
//            }
//
//            double neigh_x, neigh_y;
//            neigh_cell.get_spacial_coords(&neigh_x, &neigh_y);
//            double d = dist(x, y, neigh_x, neigh_y);
//            // cout << "dist( (" << x << "," << y << ")  (" << neigh_x << "," << neigh_y << ") )=" << d << endl;
//            double flow = - Area::a * (neigh_cell.temp[cur_temp_idx] - cell.temp[cur_temp_idx]) / d;
//            double s = get_lvl_dx(cell.lvl);
//            s = (cell.lvl <= neigh_cell.lvl) ? s : (s/2);
//
//            // проекции на оси
//            double full_flow = s * flow * fabs(x - neigh_x) / d + s * flow * fabs(y - neigh_y) / d;
//
//            termo_flows.push_back(full_flow);
//        }
//
//        double flows_sum = 0.0;
//        for (double f: termo_flows) {
//            flows_sum += f;
//        }
//
//        // вычисляем новое значение в ячейке
//
//        char border_cond_type;
//        double (*cond_func)(double, double, double);
//        stat.timers["get_border_cond"].Start();
//        cell.get_border_cond(&border_cond_type, &cond_func);
//        stat.timers["get_border_cond"].Stop();
//
//        double new_T = 0;
//
//        if (border_cond_type == char(-1)) {  // внутренняя ячейка
//
//            if (FULL_DEBUG) {
//                cout << mpiInfo.comm_rank << "cur_t=" << cell.temp[cur_temp_idx] << " flows_sum=" << flows_sum << " S=" << cell.get_S() << " Q(x,y,t)=" << Area::Q(x, y, tau * time_step_n) << endl;
//            }
//            new_T = cell.temp[cur_temp_idx] - tau * (flows_sum
//                    - Area::Q(x, y, tau * time_step_n)) /  cell.get_S();
//
//        } else if (border_cond_type == 1) {  // граничные условия первого рода
//
//            // это по сути тоже поток, только короче (половина размера ячейки) по напралению к границе
//            double border_x = x, border_y = y; // +- lvl_dx/2
//            double lvl_dx = get_lvl_dx(cell.lvl);
//            if (cell.is_right_border()) {
//                border_x += lvl_dx/2;
//            }
//            if (cell.is_left_border()) {
//                border_x -= lvl_dx/2;
//            }
//            if (cell.is_down_border()) {
//                border_y -= lvl_dx/2;
//            }
//            if (cell.is_upper_border()) {
//                border_y += lvl_dx/2;
//            }
//            double border_flow = - Area::a * (cell.temp[cur_temp_idx]
//                    - cond_func(border_x, border_y, time_step_n * tau)) / (lvl_dx/2) * lvl_dx;
//            new_T = cell.temp[cur_temp_idx] - tau * ((flows_sum + border_flow)
//                    + Area::Q(x, y, tau * time_step_n)) /  cell.get_S();
//
//        } else if (border_cond_type == 2) { // граничное условие второго рода
//
//            // TODO это просто добавить поток, заданный условием
//            double border_x = x, border_y = y; // +- lvl_dx/2
//            double lvl_dx = get_lvl_dx(cell.lvl);
//            if (cell.is_right_border()) {
//                border_x += lvl_dx/2;
//            }
//            if (cell.is_left_border()) {
//                border_x -= lvl_dx/2;
//            }
//            if (cell.is_down_border()) {
//                border_y -= lvl_dx/2;
//            }
//            if (cell.is_upper_border()) {
//                border_y += lvl_dx/2;
//            }
//            double border_flow = cond_func(border_x, border_y, time_step_n * tau);
//            new_T = cell.temp[cur_temp_idx] - tau * ((flows_sum + border_flow)
//                    + Area::Q(x, y, tau * time_step_n)) /  cell.get_S();
//
//            // std::cout << "border cond type = 2 not implemented\n";
//            // not imptemented
//        }
//
//        mesh.cells[i].temp[next_temp_idx] = new_T;
//    }

    stat.timers["step"].Stop();
}

void Proc::WriteT(string filename) {
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
                filename.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_at(fh, offset, &buf[0], len, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    stat.timers["io"].Stop();
}

void Proc::WriteStat(string filename) {
    std::cout << "MpiInfo: " << mpiInfo.toString() << std::endl;
}


void Proc::PrintMyBlocks() {

    int temp_l_corr = time_step_n % 2; // чтобы брать значение temp[0] или temp[1]
    int cur_temp_idx = temp_l_corr;

    cout << mpiInfo.comm_rank << " CELLS={ ";
    for (int i = 0; i < mesh.blocks.size(); i++) {
        cout << "(" << mesh.blocks[i].idx.lvl << "," << mesh.blocks[i].idx.i << "," << mesh.blocks[i].idx.j <<
              mesh.blocks[i].cells_lvl << ", " << mesh.blocks[i].sz << ") ";
    }
    cout << "}\n";
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
