#pragma once

#include <vector>
#include <mpi.h>
#include <math.h>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <list>

#include "../area/area.h"

using std::vector;
using std::map;
using std::string;
using std::list;

extern int    base_sz;
extern double base_dx;
extern double min_dx;
extern int    max_lvl;
extern int    base_lvl;
extern int    max_blk_lvl;
extern int    base_blk_lvl;
extern int    base_blk_sz;

void GridInit(int base_lvl, int max_lvl, int base_blk_lvl, int max_blk_lvl);

/* * * * * * * * * * РАБОТА С ИНДЕКСАМИ * * * * * * * * * */

typedef long long int GlobalNumber_t;

GlobalNumber_t merge_ints(int a, int b);
void split_ints(GlobalNumber_t c, int *a, int *b);
void split_ints(GlobalNumber_t c, int max_lvl, int *a, int *b);


enum Neigh   { DOWN, UP, LEFT, RIGHT };
enum CornerNeigh { LU, RU, LD, RD };
enum Child       { cLU, cRU, cLD, cRD };

struct TreeIndex {
    int lvl, i, j;

    TreeIndex() {}
    TreeIndex(int _lvl, int _i, int _j): lvl(_lvl), i(_i), j(_j) {}
    TreeIndex(const TreeIndex& c): lvl(c.lvl), i(c.i), j(c.j) {}
    TreeIndex(int _lvl, GlobalNumber_t globalNumber);
    TreeIndex& operator=(const TreeIndex& c) { lvl = c.lvl; i = c.i; j = c.j; return *this; }
    // TreeIndex operator=(TreeIndex c) { lvl = c.lvl; i = c.i; j = c.j; return *this; }
    // ~TreeIndex() {}

    GlobalNumber_t get_global_number();

    bool is_left_border();
    bool is_right_border();
    bool is_upper_border();
    bool is_down_border();
    bool is_left_upper_corner();
    bool is_right_upper_corner();
    bool is_left_down_corner();
    bool is_right_down_corner();

    TreeIndex get_child(Child c);
    TreeIndex get_parent();
    Child get_child_pos();

    TreeIndex get_face_neighbor(Neigh n);
    TreeIndex get_corner_neighbor(CornerNeigh n);
    vector<TreeIndex>  get_larger_possible_face_neighbour(Neigh n);
    vector<TreeIndex>& get_larger_possible_face_neighbour_optimized(vector<TreeIndex>&, Neigh n);
    vector<TreeIndex>  get_halfsize_possible_face_neighbours(Neigh n);
    vector<TreeIndex>& get_halfsize_possible_face_neighbours_optimized(vector<TreeIndex>&, Neigh n);
    
    vector<TreeIndex>  get_larger_possible_corner_neighbour(CornerNeigh n);
    vector<TreeIndex>  get_halfsize_possible_corner_neighbours(CornerNeigh n);

    vector<TreeIndex>  get_all_halfsize_possible_neighs();
    vector<TreeIndex>& get_all_halfsize_possible_neighs_optimized(vector<TreeIndex>&);
    vector<TreeIndex> get_all_samesize_possible_neighs();
    vector<TreeIndex>& get_all_samesize_possible_neighs_optimized(vector<TreeIndex>&);
    vector<TreeIndex> get_all_larger_possible_neighs();
    vector<TreeIndex>& get_all_larger_possible_neighs_optimized(vector<TreeIndex>&);

    vector<GlobalNumber_t> get_all_possible_neighbours_ids();

    bool is_border();

    void get_corner_coords(double *x, double *y);
};



/* * * * * * * * РАБОТА С ЯЧЕЙКАМИ И ДЕРЕВОМ * * * * * * * */

double dist(double x1, double y1, double x2, double y2);

struct Cell: public TreeIndex {
    char refine_mark = 0;
    double temp[2]; // for cur and next

    Cell() {}
    Cell(int _lvl, GlobalNumber_t globalNumber): TreeIndex(_lvl, globalNumber) {}
    Cell(int _lvl, int _i, int _j): TreeIndex(_lvl, _i, _j) {}
    Cell(TreeIndex ci, double _temp): TreeIndex(ci) { temp[0] = _temp; temp[1] = 0.0; }
    Cell(const Cell& c): TreeIndex(c.lvl, c.i, c.j) { temp[0] = c.temp[0]; temp[1] = c.temp[1]; refine_mark = c.refine_mark; }
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


double get_grad(double w00, double w01, double w02,
                double w10, double w11, double w12,
                double w20, double w21, double w22, double dx);
double get_lvl_dx(int lvl);





/* * * * * * * * РАБОТА БЛОЧНЫМ ДЕРЕВОМ * * * * * * * */

struct SimpleCell {
    double temp[2];
    char refine_mark = 0;

    SimpleCell() {}
    SimpleCell(int t) { temp[0] = t; temp[1] = t; }
    SimpleCell(const SimpleCell& c) { temp[0] = c.temp[0]; temp[1] = c.temp[1]; refine_mark = c.refine_mark; }
    SimpleCell& operator=(const SimpleCell& c) {
        temp[0] = c.temp[0]; temp[1] = c.temp[1]; refine_mark = c.refine_mark;
        return *this;
    }


};

struct BlockOfCells {
    vector<SimpleCell> cells;

    TreeIndex idx;
    GlobalNumber_t i;

    int cells_lvl;
    int sz;
    int refine_mark = 0;
    int refine_marks[4];

    vector<BlockOfCells*> neighs_left;
    vector<BlockOfCells*> neighs_right;
    vector<BlockOfCells*> neighs_upper;
    vector<BlockOfCells*> neighs_down;

    vector<GlobalNumber_t> neighs_left_idxs;
    vector<GlobalNumber_t> neighs_right_idxs;
    vector<GlobalNumber_t> neighs_upper_idxs;
    vector<GlobalNumber_t> neighs_down_idxs;

    size_t GetMemSize() {
        size_t res = 0;
        res += sizeof(SimpleCell) * cells.size();
        res += sizeof(BlockOfCells*) * neighs_left.size();
        res += sizeof(BlockOfCells*) * neighs_right.size();
        res += sizeof(BlockOfCells*) * neighs_upper.size();
        res += sizeof(BlockOfCells*) * neighs_down.size();
        res += sizeof(GlobalNumber_t) * neighs_left_idxs.size();
        res += sizeof(GlobalNumber_t) * neighs_right_idxs.size();
        res += sizeof(GlobalNumber_t) * neighs_upper_idxs.size();
        res += sizeof(GlobalNumber_t) * neighs_down_idxs.size();

        return res;
    }

    BlockOfCells(int _cells_lvl, int _blk_lvl, int _sz, GlobalNumber_t _i);
    BlockOfCells(TreeIndex _idx, int _cells_lvl, int _sz): idx(_idx), i(_idx.get_global_number()), cells_lvl(_cells_lvl), sz(_sz) {
        refine_marks[0] = 0;
        refine_marks[1] = 0;
        refine_marks[2] = 0;
        refine_marks[3] = 0;
    }
    BlockOfCells(const BlockOfCells& b) {
        idx = b.idx; i = b.i;
        cells = b.cells;
        cells_lvl = b.cells_lvl; sz = b.sz; refine_mark = b.refine_mark; 
        refine_marks[0] = b.refine_marks[0]; refine_marks[1] = b.refine_marks[1]; 
        refine_marks[2] = b.refine_marks[2]; refine_marks[3] = b.refine_marks[3];
        neighs_left_idxs = b.neighs_left_idxs;
        neighs_right_idxs = b.neighs_right_idxs;
        neighs_upper_idxs = b.neighs_upper_idxs;
        neighs_down_idxs = b.neighs_down_idxs;
    }
    BlockOfCells& operator=(const BlockOfCells& b) {
        idx = b.idx; i = b.i;
        cells = b.cells;
        cells_lvl = b.cells_lvl; sz = b.sz; refine_mark = b.refine_mark; 
        refine_marks[0] = b.refine_marks[0]; refine_marks[1] = b.refine_marks[1]; 
        refine_marks[2] = b.refine_marks[2]; refine_marks[3] = b.refine_marks[3];
        return *this;
    }

    vector<BlockOfCells> Split(double (*Temp_func)(double, double));

    int MarkToRefine();
    void mark_quarter(int i, int j);
    void RefineCells();
    void ClearMarks();

    // реально создать ячейки
    void CreateCells(double (*Temp_func)(double, double));

    int GetNOfCells() { return cells.size(); }

    void get_spacial_coords(int i, int j, double *x, double *y);
    SimpleCell *find_border_cell_by_global_idx(GlobalNumber_t global_idx);

};

GlobalNumber_t get_glob_idx(GlobalNumber_t blk_i, GlobalNumber_t cell_i, int cell_lvl);
int check_cell_owner(GlobalNumber_t blk_i, GlobalNumber_t cell_i);
vector<GlobalNumber_t> find_cell_neighs_ids_in_blk(GlobalNumber_t cell_glob_idx, int cell_lvl, BlockOfCells* blk, Neigh neigh_dir);

struct BlockedLinearTree {
    int max_present_lvl;
    int max_present_blk_lvl;
    // vector<int> block_offsets;
    vector<BlockOfCells> blocks;  // это линейное дерево
    vector<int> proc_blocks;


    BlockedLinearTree() { max_present_lvl = 0; max_present_blk_lvl = 0; }
    BlockedLinearTree(double (*Temp_func)(double, double));


    // returns 1 if there are cells/blocks to refine, 0 otherwise
    int MarkToRefine();

    void RefineBlocks();
    void RefineCells();

    void BuildNeighs();
    BlockOfCells* find_block(TreeIndex target);
    BlockOfCells* find_block(GlobalNumber_t target);
    vector<BlockOfCells*> find_block_children(TreeIndex target, Neigh n);

    void Decompose(int n_procs);

    // Write for visualization
    void Write(string filename);

    // TODO i/o
    void WriteBlocks(string filename);
    void WriteOffsets(string filename);
    void GenFromWriteBlocksStruct(vector<char> buf, double (*start_func)(double, double));
//    void ReadBlocks();


    vector<char> GenWriteStruct(char light);

private:
    BlockedLinearTree(const BlockedLinearTree&);
    BlockedLinearTree& operator=(const BlockedLinearTree&);


    vector<char> GenWriteBlocksStruct();
};

struct WriteCell {
    int lvl;
    int i, j;
    double temp;
    int proc;
    int blk;
};

struct WriteBlock {
    int lvl, i, j;
    int cells_lvl;
    int sz;
};
