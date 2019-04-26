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

    if (ghosts_out_ids != nullptr) {
        delete[] ghosts_out_ids;
    }
    if (ghosts_in != nullptr) {
        delete[] ghosts_in;
    }
    if (ghosts_in_temps != nullptr) {
        delete[] ghosts_in_temps;
    }
    if (ghosts_out_temps != nullptr) {
        delete[] ghosts_out_temps;
    }

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
    stat.timers["build_ghosts"].Start();

    MPI_Barrier(mpiInfo.comm);

    ghosts_in_temps  = new vector<double>[mpiInfo.comm_size];
    ghosts_out_temps = new vector<double>[mpiInfo.comm_size];

    ghosts_out_ids = new vector<GlobalNumber_t>[mpiInfo.comm_size];
    ghosts_in      = new LinearTree[mpiInfo.comm_size];

    // cout << mpiInfo.comm_rank << " BuildGhosts here1\n";

    for (Cell c: mesh.cells) {
        vector<Cell> neighs;
        
        vector<GlobalNumber_t> neigh_ids = c.get_all_possible_neighbours_ids();

        if (FULL_DEBUG) {
            cout << mpiInfo.comm_rank << " neigh_ids(" << c.get_global_number() << ")= {sz=" << neigh_ids.size() << "}";
            for (int i = 0; i < neigh_ids.size(); i++) {
                cout << neigh_ids[i] << " ";
            }
            cout << endl;
        }

        // find their owners
        for (GlobalNumber_t neigh: neigh_ids) {
            int owner = find_owner(neigh);
            if (owner == -1) {
                continue;
            }
            if (owner != mpiInfo.comm_rank) {
                ghosts_out_ids[owner].push_back(c.get_global_number());
            }
        }
    }

    // cout << mpiInfo.comm_rank << " BuildGhosts here2\n";

    // сортируем и удаляем повторения
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        std::sort(ghosts_out_ids[i].begin(), ghosts_out_ids[i].end());
        auto last = std::unique(ghosts_out_ids[i].begin(), ghosts_out_ids[i].end());
        ghosts_out_ids[i].erase(last, ghosts_out_ids[i].end());
    }

    if (FULL_DEBUG) {
        cout << mpiInfo.comm_rank << " GHOSTS_OUT_IDS={ ";
        for (int n = 0; n < mpiInfo.comm_size; n++) {
            cout << "[ ";
            for (int i = 0; i < ghosts_out_ids[n].size(); i++) {
                cout << ghosts_out_ids[n][i] << " ";
            }
            cout << "] ";
        }
        cout << "}\n";
    }

    // cout << mpiInfo.comm_rank << " BuildGhosts here3\n";

    // считаем количество соседей, которым что-то нужно отправлять и обмениваемся размерами
    vector<int> out_lens(mpiInfo.comm_size, 0);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        out_lens[i] = ghosts_out_ids[i].size();
    }
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        if (out_lens[i] > 0) {
            active_neighs_num++;
        }
    }
    cout << mpiInfo.comm_rank << " active_neighs_num=" << active_neighs_num << endl;
    vector<int> in_lens(mpiInfo.comm_size, 0);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        MPI_Gather(&out_lens[i], 1, MPI_INT, &in_lens[0], 1, MPI_INT, i, mpiInfo.comm);
    }

    if (FULL_DEBUG) {
        cout << mpiInfo.comm_rank << " IN_LENS=[ ";
        for (int i = 0; i < mpiInfo.comm_size; i++) {
            cout << in_lens[i] << " ";
        }
        cout << "]\n";
    }
    
    // cout << mpiInfo.comm_rank << " BuildGhosts here4\n";

    // обмениваемся idшниками ячеек от соседей и формируем LinearTree ячеек
    vector<vector<GlobalNumber_t>> ghosts_in_ids_tmp(mpiInfo.comm_size);
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        ghosts_in_ids_tmp[i] = vector<GlobalNumber_t>(in_lens[i]);
    }

    // cout << mpiInfo.comm_rank << " BuildGhosts here44\n";

    MPI_Request send_reqs[active_neighs_num];
    MPI_Status send_statuses[active_neighs_num];
    int req_num = 0;

    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if (out_lens[n] > 0) {
            MPI_Isend(&ghosts_out_ids[n][0], out_lens[n], MPI_LONG_LONG_INT, n, 0, mpiInfo.comm, &send_reqs[req_num]);
            req_num++;
        }
    }

    // cout << mpiInfo.comm_rank << " BuildGhosts here444\n";

    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ( (n != mpiInfo.comm_rank) && (in_lens[n] > 0) ) {
            MPI_Status status;
            MPI_Recv(&ghosts_in_ids_tmp[n][0], in_lens[n], MPI_LONG_LONG_INT, n, MPI_ANY_TAG, mpiInfo.comm, &status);

            // fill temp in in_ghosts[n]
            for (int i = 0; i < in_lens[n]; i++) {
                ghosts_in[n].cells.push_back(Cell(0, ghosts_in_ids_tmp[n][i])); 
                ghosts_in[n].cells[i].temp[0] = 0.0;
            }
        }
    }

    // cout << mpiInfo.comm_rank << " BuildGhosts here4444\n";

    MPI_Waitall(active_neighs_num, send_reqs, send_statuses);

    if (FULL_DEBUG) {
        cout << mpiInfo.comm_rank << " GHOSTS_IN_WITHOUT_LEVELS={ ";
        for (int n = 0; n < mpiInfo.comm_size; n++) {
            cout << "[ ";
            for (int i = 0; i < ghosts_in[n].cells.size(); i++) {
                cout << "(" << ghosts_in[n].cells[i].lvl << "," << ghosts_in[n].cells[i].i << "," << ghosts_in[n].cells[i].j << ") ";
            }
            cout << "] ";
        }
        cout << "}\n";
    }

    // cout << mpiInfo.comm_rank << " BuildGhosts here5\n";

    // заполняем уровни ячеек
    vector<vector<int>> ghosts_out_lvls_tmp(mpiInfo.comm_size);
    for (int n = 0; n < mpiInfo.comm_size; n++ ) {
        for (int i = 0; i < out_lens[n]; i++) {
            Cell c;
            int c_idx = mesh.FindCell(ghosts_out_ids[n][i], &c);
            if (c_idx == -1) {
                cout << "ERROR no cell found\n";
                // continue;
            }
            int c_lvl = mesh.cells[c_idx].lvl;
            ghosts_out_lvls_tmp[n].push_back(c_lvl);
        }
    }

    if (FULL_DEBUG) {
        cout << mpiInfo.comm_rank << " GHOSTS_OUT_LEVELS={ ";
        for (int n = 0; n < mpiInfo.comm_size; n++) {
            cout << "[ ";
            for (int i = 0; i < ghosts_out_lvls_tmp[n].size(); i++) {
                cout << ghosts_out_lvls_tmp[n][i] << " ";
            }
            cout << "] ";
        }
        cout << "}\n";
    }

    req_num = 0;
    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if (out_lens[n] > 0) {
            MPI_Isend(&ghosts_out_lvls_tmp[n][0], out_lens[n], MPI_INT, n, 0, mpiInfo.comm, &send_reqs[req_num]);
            req_num++;
        }
    }

    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ( (n != mpiInfo.comm_rank) && (in_lens[n] > 0) ) {
            MPI_Status status;
            vector<int>buf(in_lens[n]);
            MPI_Recv(&buf[0], in_lens[n], MPI_INT, n, MPI_ANY_TAG, mpiInfo.comm, &status);

            // fill lvlvs in in_ghosts[n]
            for (int i = 0; i < ghosts_in[n].cells.size(); i++) {
                ghosts_in[n].cells[i].lvl = buf[i];
            }
        }
    }

    MPI_Waitall(active_neighs_num, send_reqs, send_statuses);

    // cout << mpiInfo.comm_rank << " BuildGhosts here6\n";

    // создаем буферы для входящих и исходящих температур
    for (int i = 0; i < mpiInfo.comm_size; i++) {
        ghosts_in_temps[i]  = vector<double>(in_lens[i], 0);
        ghosts_out_temps[i] = vector<double>(out_lens[i], 0);
    }

    // cout << mpiInfo.comm_rank << " BuildGhosts here7\n";

    std::cout << mpiInfo.comm_rank << " ghosts buffers built\n";

    if (FULL_DEBUG) {
        cout << mpiInfo.comm_rank << " GHOSTS={ ";
        for (int n = 0; n < mpiInfo.comm_size; n++) {
            cout << "[ ";
            for (int i = 0; i < ghosts_in[n].cells.size(); i++) {
                cout << "(" << ghosts_in[n].cells[i].lvl << "," << ghosts_in[n].cells[i].i << "," << ghosts_in[n].cells[i].j << ") ";
            }
            cout << "] ";
        }
        cout << "}\n";
    }


    stat.timers["build_ghosts"].Stop();
    return 0;
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

int Proc::ExchangeGhosts() {
    
    // cout << mpiInfo.comm_rank << " ExchangeGhosts started\n";
    stat.timers["exchange_ghosts"].Start();

    int temp_l_corr = time_step_n % 2; // чтобы брать значение temp[0] или temp[1]
    int cur_temp_idx = temp_l_corr;

    // non-blocking send to all neighs
    MPI_Request send_reqs[active_neighs_num];
    MPI_Status send_statuses[active_neighs_num];
    int req_num = 0;

    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ( (n != mpiInfo.comm_rank) && (ghosts_out_ids[n].size() > 0) ) {
            // fill out temps for n
            for (int i = 0; i < ghosts_out_ids[n].size(); i++) {
                Cell c;
                int idx = mesh.FindCell(ghosts_out_ids[n][i], &c);
                if (idx == -1) {
                    cout << "[ERROR] exchange ghosts cell not found but have to\n";
                }
                ghosts_out_temps[n][i] = mesh.cells[idx].temp[cur_temp_idx];
            }
            
            MPI_Isend(&ghosts_out_temps[n][0], ghosts_out_ids[n].size(), MPI_DOUBLE, n, 0, mpiInfo.comm, &send_reqs[req_num]);
            req_num++;
        }
    }

    // blocking recv from all neighs

    for (int n = 0; n < mpiInfo.comm_size; n++) {
        if ( (n != mpiInfo.comm_rank) && (ghosts_in[n].cells.size() > 0) ) {
            MPI_Status status;
            MPI_Recv(&ghosts_in_temps[n][0], ghosts_in[n].cells.size(), MPI_DOUBLE, n, MPI_ANY_TAG, mpiInfo.comm, &status);

            // fill temp in in_ghosts[n]
            for (int i = 0; i < ghosts_in[n].cells.size(); i++) {
                ghosts_in[n].cells[i].temp[cur_temp_idx] = ghosts_in_temps[n][i];
            }
        }
    }

    
    MPI_Waitall(active_neighs_num, send_reqs, send_statuses);


    if (FULL_DEBUG) {
        cout << mpiInfo.comm_rank << " GHOSTS={ ";
        for (int n = 0; n < mpiInfo.comm_size; n++) {
            cout << "[ ";
            for (int i = 0; i < ghosts_in[n].cells.size(); i++) {
                cout << "(" << ghosts_in[n].cells[i].lvl << "," << ghosts_in[n].cells[i].i << "," << ghosts_in[n].cells[i].j << ": " <<  ghosts_in[n].cells[i].temp[cur_temp_idx] << ") ";
            }
            cout << "] ";
        }
        cout << "}\n";
    }


    stat.timers["exchange_ghosts"].Stop();
    // cout << mpiInfo.comm_rank << " ExchangeGhosts finished\n";
    return 0;
}

void Proc::MakeStep() {

    // usleep(10000000);

    stat.timers["step"].Start();
    int temp_l_corr = time_step_n % 2; // чтобы брать значение temp[0] или temp[1]
    int cur_temp_idx = temp_l_corr;
    int next_temp_idx = (temp_l_corr + 1) % 2;

    ExchangeGhosts();

    // PrintMyCells();
    // PrintGhostCells();

    time_step_n++;

    for (int i = 0; i < mesh.cells.size(); i++) {
        
        Cell cell = mesh.cells[i];
        
        double x, y;
        cell.get_spacial_coords(&x, &y);

        // ищем всех возможных соседей
        stat.timers["get_possible_neighs"].Start();
        vector<GlobalNumber_t> possible_neigh_ids = cell.get_all_possible_neighbours_ids();
        stat.timers["get_possible_neighs"].Stop();


        // вычисляем потоки междуу ячейкой и реально существующими соседями
        vector<double> termo_flows;
        for (GlobalNumber_t id: possible_neigh_ids) {
            int owner = find_owner(id);
            if (owner == -1) {
                continue;
            }
            Cell neigh_cell;

            if (owner == mpiInfo.comm_rank) { // ячейка у меня
                stat.timers["find_cell"].Start();
                int idx = mesh.FindCell(id, &neigh_cell);
                stat.timers["find_cell"].Stop();
                if (idx == -1) {
                    continue;
                }
            } else { // ячейка-призрак
                stat.timers["find_cell"].Start();
                int idx = ghosts_in[owner].FindCell(id, &neigh_cell);
                stat.timers["find_cell"].Stop();
                if (idx == -1) {
                    continue;
                }
            }

            double neigh_x, neigh_y;
            neigh_cell.get_spacial_coords(&neigh_x, &neigh_y);
            double d = dist(x, y, neigh_x, neigh_y);
            // cout << "dist( (" << x << "," << y << ")  (" << neigh_x << "," << neigh_y << ") )=" << d << endl;
            double flow = - Area::a * (neigh_cell.temp[cur_temp_idx] - cell.temp[cur_temp_idx]) / d;
            double s = get_lvl_dx(cell.lvl);
            s = (cell.lvl <= neigh_cell.lvl) ? s : (s/2);

            // проекции на оси
            double full_flow = s * flow * fabs(x - neigh_x) / d + s * flow * fabs(y - neigh_y) / d;

            termo_flows.push_back(full_flow);
        }

        double flows_sum = 0.0;
        for (double f: termo_flows) {
            flows_sum += f;
        }

        // вычисляем новое значение в ячейке
       
        char border_cond_type;
        double (*cond_func)(double, double, double);
        stat.timers["get_border_cond"].Start();
        cell.get_border_cond(&border_cond_type, &cond_func);
        stat.timers["get_border_cond"].Stop();

        double new_T = 0;
        
        if (border_cond_type == char(-1)) {  // внутренняя ячейка

            if (FULL_DEBUG) {
                cout << mpiInfo.comm_rank << "cur_t=" << cell.temp[cur_temp_idx] << " flows_sum=" << flows_sum << " S=" << cell.get_S() << " Q(x,y,t)=" << Area::Q(x, y, tau * time_step_n) << endl;
            }
            new_T = cell.temp[cur_temp_idx] - tau * (flows_sum  
                    - Area::Q(x, y, tau * time_step_n)) /  cell.get_S();

        } else if (border_cond_type == 1) {  // граничные условия первого рода

            // это по сути тоже поток, только короче (половина размера ячейки) по напралению к границе
            double border_x = x, border_y = y; // +- lvl_dx/2
            double lvl_dx = get_lvl_dx(cell.lvl);
            if (cell.is_right_border()) {
                border_x += lvl_dx/2;
            }
            if (cell.is_left_border()) {
                border_x -= lvl_dx/2;
            }
            if (cell.is_down_border()) {
                border_y -= lvl_dx/2;
            }
            if (cell.is_upper_border()) {
                border_y += lvl_dx/2;
            }
            double border_flow = - Area::a * (cell.temp[cur_temp_idx] 
                    - cond_func(border_x, border_y, time_step_n * tau)) / (lvl_dx/2) * lvl_dx;
            new_T = cell.temp[cur_temp_idx] - tau * ((flows_sum + border_flow) 
                    + Area::Q(x, y, tau * time_step_n)) /  cell.get_S();

        } else if (border_cond_type == 2) { // граничное условие второго рода

            // TODO это просто добавить поток, заданный условием
            double border_x = x, border_y = y; // +- lvl_dx/2
            double lvl_dx = get_lvl_dx(cell.lvl);
            if (cell.is_right_border()) {
                border_x += lvl_dx/2;
            }
            if (cell.is_left_border()) {
                border_x -= lvl_dx/2;
            }
            if (cell.is_down_border()) {
                border_y -= lvl_dx/2;
            }
            if (cell.is_upper_border()) {
                border_y += lvl_dx/2;
            }
            double border_flow = cond_func(border_x, border_y, time_step_n * tau);
            new_T = cell.temp[cur_temp_idx] - tau * ((flows_sum + border_flow) 
                    + Area::Q(x, y, tau * time_step_n)) /  cell.get_S();

            // std::cout << "border cond type = 2 not implemented\n";
            // not imptemented
        }

        mesh.cells[i].temp[next_temp_idx] = new_T;
    }
    stat.timers["step"].Stop();
}

void Proc::WriteT(string filename) {
    stat.timers["io"].Start();

    vector<char> buf = mesh.GenWriteStruct();
    int len = buf.size();
    int *lens = new int[mpiInfo.comm_size];
    stat.timers["communication"].Start();
    MPI_Allgather(&len, 1, MPI_INT, (void *)lens, 1, MPI_INT, mpiInfo.comm);
    stat.timers["communication"].Stop();
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

int Proc::ISendGhosts() {
    // MPI_ISend()

    return 0;
}
int Proc::IRecvGhosts() {
    // MPI_IRecv()
    // запомнить шутку, по которой ждать потом
    return 0;
}
int Proc::WaitallGhosts() {
    return 0;
}

// void MPI_Allgather_wrapper()

void Proc::PrintMyCells() {

    int temp_l_corr = time_step_n % 2; // чтобы брать значение temp[0] или temp[1]
    int cur_temp_idx = temp_l_corr;

    cout << mpiInfo.comm_rank << " CELLS={ ";
    for (int i = 0; i < mesh.cells.size(); i++) {
        cout << "(" << mesh.cells[i].lvl << "," << mesh.cells[i].i << "," << mesh.cells[i].j << ": " <<  mesh.cells[i].temp[cur_temp_idx] << ") ";
    }
    cout << "}\n";
}

void Proc::PrintGhostCells() {

    int temp_l_corr = time_step_n % 2; // чтобы брать значение temp[0] или temp[1]
    int cur_temp_idx = temp_l_corr;

    cout << mpiInfo.comm_rank << " GHOSTS={ ";
        for (int n = 0; n < mpiInfo.comm_size; n++) {
            cout << n << ":[ ";
            for (int i = 0; i < ghosts_in[n].cells.size(); i++) {
                cout << "(" << ghosts_in[n].cells[i].lvl << "," << ghosts_in[n].cells[i].i << "," << ghosts_in[n].cells[i].j << ": " <<  ghosts_in[n].cells[i].temp[cur_temp_idx] << ") ";
            }
            cout << "] ";
        }
        cout << "}\n";
}
