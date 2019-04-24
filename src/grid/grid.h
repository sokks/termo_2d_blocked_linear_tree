#pragma once

#include <vector>
#include <mpi.h>
#include <math.h>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>

#include "../area/area.h"

using std::vector;
using std::map;
using std::string;

extern int    base_sz;
extern double base_dx;
extern double min_dx;
extern int    max_lvl;
extern int    base_lvl;

void GridInit(int base_lvl, int max_lvl, int base_blk_lvl, int max_blk_lvl);

/* * * * * * * * * * РАБОТА С ИНДЕКСАМИ * * * * * * * * * */

typedef long long int GlobalNumber_t;

GlobalNumber_t merge_ints(int a, int b);
void split_ints(GlobalNumber_t c, int *a, int *b);


enum Neigh   { DOWN, UP, LEFT, RIGHT };
enum CornerNeigh { LU, RU, LD, RD };
enum Child       { cLU, cRU, cLD, cRD };

struct CellIndex {
    int lvl, i, j;

    CellIndex() {}
    CellIndex(int _lvl, int _i, int _j): lvl(_lvl), i(_i), j(_j) {}
    CellIndex(const CellIndex& c): lvl(c.lvl), i(c.i), j(c.j) {}
    CellIndex(int _lvl, GlobalNumber_t globalNumber);
    CellIndex& operator=(const CellIndex& c) { lvl = c.lvl; i = c.i; j = c.j; return *this; }
    // CellIndex operator=(CellIndex c) { lvl = c.lvl; i = c.i; j = c.j; return *this; }
    // ~CellIndex() {}

    GlobalNumber_t get_global_number();

    bool is_left_border();
    bool is_right_border();
    bool is_upper_border();
    bool is_down_border();
    bool is_left_upper_corner();
    bool is_right_upper_corner();
    bool is_left_down_corner();
    bool is_right_down_corner();

    CellIndex get_child(Child c);
    CellIndex get_parent();
    Child get_child_pos();

    CellIndex get_face_neighbor(Neigh n);
    CellIndex get_corner_neighbor(CornerNeigh n);
    vector<CellIndex>  get_larger_possible_face_neighbour(Neigh n);
    vector<CellIndex>& get_larger_possible_face_neighbour_optimized(vector<CellIndex>&, Neigh n);
    vector<CellIndex>  get_halfsize_possible_face_neighbours(Neigh n);
    vector<CellIndex>& get_halfsize_possible_face_neighbours_optimized(vector<CellIndex>&, Neigh n);
    
    vector<CellIndex>  get_larger_possible_corner_neighbour(CornerNeigh n);
    vector<CellIndex>  get_halfsize_possible_corner_neighbours(CornerNeigh n);

    vector<CellIndex>  get_all_halfsize_possible_neighs();
    vector<CellIndex>& get_all_halfsize_possible_neighs_optimized(vector<CellIndex>&);
    vector<CellIndex> get_all_samesize_possible_neighs();
    vector<CellIndex>& get_all_samesize_possible_neighs_optimized(vector<CellIndex>&);
    vector<CellIndex> get_all_larger_possible_neighs();
    vector<CellIndex>& get_all_larger_possible_neighs_optimized(vector<CellIndex>&);

    vector<GlobalNumber_t> get_all_possible_neighbours_ids();

    bool is_border();
};



/* * * * * * * * РАБОТА С ЯЧЕЙКАМИ И ДЕРЕВОМ * * * * * * * */

double dist(double x1, double y1, double x2, double y2);

struct Cell: public CellIndex {
    char refine_mark = 0;
    double temp[2]; // for cur and next

    Cell() {}
    Cell(int _lvl, GlobalNumber_t globalNumber): CellIndex(_lvl, globalNumber) {}
    Cell(int _lvl, int _i, int _j): CellIndex(_lvl, _i, _j) {}
    Cell(CellIndex ci, double _temp): CellIndex(ci) { temp[0] = _temp; temp[1] = 0.0; }
    Cell(const Cell& c): CellIndex(c.lvl, c.i, c.j) { temp[0] = c.temp[0]; temp[1] = c.temp[1]; refine_mark = c.refine_mark; }
    Cell& operator=(const Cell& c) { lvl = c.lvl; i = c.i; j = c.j; temp[0] = c.temp[0]; temp[1] = c.temp[1]; refine_mark = c.refine_mark; return *this; }

    double *get_temp() { return &temp[0]; }
    void    get_spacial_coords(double *x, double *y);
    void    get_border_cond(char *cond_type, double (**cond_func)(double, double, double));
    double  get_S() {
        double lvl_dx = min_dx * pow(2, max_lvl - lvl);
        return lvl_dx * lvl_dx;
    }

    void mark_to_refine() { refine_mark = 1; }
    vector<Cell> split();

};


// ! очень большая оптимизация по времени и памяти, если искать в диапазоне

struct LinearTree {
    
    int max_present_lvl;
    vector<Cell> cells;

    LinearTree(int reserved_size=0) { cells = vector<Cell>(reserved_size); max_present_lvl = 0; }
    LinearTree(string filename);
    LinearTree(double (*Temp_func)(double, double));

    int FindCell(GlobalNumber_t target, Cell *cell);

    // MarkToRefine returns 1 if where are cells to refine and they can be refined, 0 otherwise
    int  MarkToRefine();
    void DoRefine();
    void Balance21();

    void Write(string filename);
    void WriteOffsets(string filename, int n_of_procs);
    vector<char> GenWriteStruct();
    void GenFromWriteStruct(vector<char>&);

};

double get_grad(double w00, double w01, double w02,
                double w10, double w11, double w12,
                double w20, double w21, double w22, double dx);
double get_lvl_dx(int lvl);





/* * * * * * * * РАБОТА БЛОЧНЫМ ДЕРЕВОМ * * * * * * * */

struct BlockOfCells {
    vector<Cell> cells;
    GlobalNumber_t i;
    int cells_lvl;
    int blk_lvl;
    int sz;
    int refine_mark = 0;

    
    BlockOfCells(int _cells_lvl, int _blk_lvl, int _sz, GlobalNumber_t _i): 
        cells_lvl(_cells_lvl), _blk_lvl(blk_lvl), sz(_sz), i(_i) { 
        // blk = vector<Cell>(_sz); 
    }

    // TODO реально создать ячейки
    void CreateCells(double (*Temp_func)(double, double));

    int GetSize() { return blk.size(); }

    // TODO пройти по ячейкам и если хотя бы одна из них должна быть измельчена, пометить блок
    void AnalyseRefine() { }
};


struct BlockedLinearTree {
    int block_size = 0;
    int max_present_lvl;
    int max_blk_lvl;
    // vector<int> block_offsets;
    vector<BlockOfCells> blocks;



    BlockedLinearTree(int reserved_size=0) { blocks = vector<BlockOfCells>(reserved_size); max_present_lvl = 0; max_blk_lvl = 0; }
    BlockedLinearTree(double (*Temp_func)(double, double));
 
    // строим блоки: покрываем разбиение, записываем оффсеты
    // возвращает процент покрытия
    double BuildBlocks(int block_sz = 16);

    // перегружаем на более простой и быстрый вариант на основе блоков
    // не должно быть обращений к памяти ячеек
    int FindCell(GlobalNumber_t target, Cell *cell);

    // здесь изменяется стратегия разбиения, можно только по концам блоков
    void WriteOffsets(string filename, int n_of_procs);
};
