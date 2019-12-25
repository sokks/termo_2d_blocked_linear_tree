#include "grid.h"

using std::vector;
using std::map;
using std::string;
using std::sort;
using std::cout;
using std::endl;

int    base_sz;
double base_dx;
double min_dx;
int    max_lvl;
int    base_lvl;
int    base_blk_lvl;
int    max_blk_lvl;
int    base_blk_sz; // ширина = высота в ячейках (количество ячеек в блоке = base_blk_sz^2)

void GridInit(int base_level, int max_level, int base_blk_level, int max_blk_level) {
    base_lvl = base_level;
    max_lvl  = max_level;
    base_blk_lvl = base_blk_level;
    max_blk_lvl  = max_blk_level;
    
    base_sz  = pow(2, base_lvl);
    base_dx  = (Area::x_end - Area::x_start) / base_sz;
    int n_blocks_row = pow(2, base_blk_lvl);
    base_blk_sz  = base_sz / n_blocks_row;

    int max_sz = pow(2, max_lvl);
    min_dx  = (Area::x_end - Area::x_start) / max_sz;

    cout << "GRID PARAMS INITED" <<
            "\nbase_lvl=" << base_lvl <<
            "\nmax_lvl=" << max_lvl <<
            "\nbase_blk_lvl=" << base_blk_lvl <<
            "\nmax_blk_lvl=" << max_blk_lvl <<
            "\nbase_sz=" << base_sz << endl <<
            "\nbase_dx=" << base_dx <<
            "\nbase_blk_sz=" << base_blk_sz <<
            "\nmin_dx=" << min_dx << endl;
}


/* * * * * * * * * * РАБОТА С ИНДЕКСАМИ * * * * * * * * * */

GlobalNumber_t merge_ints(int a, int b) {
    GlobalNumber_t c = 0;

    int c_pos = max_lvl * 2;
    int a_pos = max_lvl;
    int b_pos = max_lvl;

    for (int pos = max_lvl - 1; pos >= 0; pos--) {
        int a1 = a & (1 << pos);
        GlobalNumber_t a2 = a1 << (pos + 1);
        c = c | a2;
        
        int b1 = b & (1 << pos);
        GlobalNumber_t b2 = b1 << pos;
        c = c | b2;
    }

    return c;
}

void split_ints(GlobalNumber_t c, int *a, int *b) {
    int aa = 0;
    int bb = 0;

    for (int pos = max_lvl - 1; pos >= 0; pos--) {
        GlobalNumber_t c1 = c & (1 << (pos*2+1));
        c1 = c1 >> (pos + 1);
        aa = aa | c1;

        GlobalNumber_t c2 = c & (1 << (pos*2));
        c2 = c2 >> (pos);
        bb = bb | c2;
    }

    *a = aa;
    *b = bb;
}

void split_ints(GlobalNumber_t c, int max_lvl, int *a, int *b) {
    int aa = 0;
    int bb = 0;

    for (int pos = max_lvl - 1; pos >= 0; pos--) {
        GlobalNumber_t c1 = c & (1 << (pos*2+1));
        c1 = c1 >> (pos + 1);
        aa = aa | c1;

        GlobalNumber_t c2 = c & (1 << (pos*2));
        c2 = c2 >> (pos);
        bb = bb | c2;
    }

    *a = aa;
    *b = bb;
}

TreeIndex::TreeIndex(int _lvl, GlobalNumber_t globalNumber): lvl(_lvl) {
    // std::cout << "new TreeIndex(" << lvl << ", " << globalNumber << ")\n";
    split_ints(globalNumber, &i, &j);
}

GlobalNumber_t TreeIndex::get_global_number() {
    // cout << "global_number(" << lvl << "," << i << "," << j << ")=" << merge_ints(i, j) << endl;
    return merge_ints(i, j);
}

TreeIndex TreeIndex::get_child(Child c) {
    TreeIndex child;
    child.lvl = lvl+1;
    int h = 1 << (max_lvl - child.lvl);

    if (c == cLD) {        // (0,0)
        child.i = i;
        child.j = j;
    } else if (c == cRD) { // (0,1)
        child.i = i;
        child.j = j + h;
    } else if (c == cLU) { // (1,0)
        child.i = i + h;
        child.j = j;
    } else {                     // (1,1)
        child.i = i + h;
        child.j = j + h;
    }

    return child;
}

TreeIndex TreeIndex::get_parent() {
    TreeIndex parent;
    parent.lvl = lvl - 1;
    parent.i = i & ~(1 << (max_lvl - lvl));
    parent.j = j & ~(1 << (max_lvl - lvl));
    return parent;
}

Child TreeIndex::get_child_pos() {
    int hi = i & (1 << (max_lvl - lvl));
    int hj = j & (1 << (max_lvl - lvl));
    if ((hi == 0) && (hj == 0)) {
        return cLD;
    }
    if ((hi == 0) && (hj != 0)) {
        return cRD;
    }
    if ((hi != 0) && (hj == 0)) {
        return cLU;
    }
    if ((hi != 0) && (hj != 0)) {
        return cRU;
    }
    return cLD;
}

TreeIndex TreeIndex::get_face_neighbor(Neigh n) {
    int h = 1 << (max_lvl - lvl); // 2^(b-l)
    TreeIndex c;
    c.lvl = lvl;
    c.i   = i + ((n == DOWN) ? -h : (n == UP) ? h : 0);
    c.j   = j + ((n == LEFT) ? -h : (n == RIGHT) ? h : 0);
    return c;
}

TreeIndex TreeIndex::get_corner_neighbor(CornerNeigh n) {
    int h = 1 << (max_lvl - lvl - 1); // 2^(b-l)
    TreeIndex c;
    c.lvl = lvl;
    c.i   = i + ( ((n == LU) || (n == LD)) ? -h : h);
    c.j   = j + ( ((n == LD) || (n == RD)) ? -h : h);
    return c;
}

vector<TreeIndex> TreeIndex::get_larger_possible_face_neighbour(Neigh n) {
    vector<TreeIndex> ret;

    Child my_pos = get_child_pos();
    if (n == RIGHT) {
        if ((my_pos == cLD) || (my_pos == cLU)) {
            return ret;
        }
        TreeIndex c = get_parent();
        ret.push_back(c.get_face_neighbor(n));
        return ret;
    }

    if (n == LEFT) {
        if ((my_pos == cRD) || (my_pos == cRU)) {
            return ret;
        }
        TreeIndex c = get_parent();
        ret.push_back(c.get_face_neighbor(n));
        return ret;
    }

    if (n == UP) {
        if ((my_pos == cLD) || (my_pos == cRD)) {
            return ret;
        }
        TreeIndex c = get_parent();
        ret.push_back(c.get_face_neighbor(n));
        return ret;
    }
    
    if (n == DOWN) {
        if ((my_pos == cLU) || (my_pos == cRU)) {
            return ret;
        }
        TreeIndex c = get_parent();
        ret.push_back(c.get_face_neighbor(n));
        return ret;
    }

    return ret;
}

vector<TreeIndex>& TreeIndex::get_larger_possible_face_neighbour_optimized(vector<TreeIndex>& buf, Neigh n) {
    Child my_pos = get_child_pos();
    if (n == RIGHT) {
        if ((my_pos == cLD) || (my_pos == cLU)) {
            return buf;
        }
        TreeIndex c = get_parent();
        buf.push_back(c.get_face_neighbor(n));
        return buf;
    }

    if (n == LEFT) {
        if ((my_pos == cRD) || (my_pos == cRU)) {
            return buf;
        }
        TreeIndex c = get_parent();
        buf.push_back(c.get_face_neighbor(n));
        return buf;
    }

    if (n == UP) {
        if ((my_pos == cLD) || (my_pos == cRD)) {
            return buf;
        }
        TreeIndex c = get_parent();
        buf.push_back(c.get_face_neighbor(n));
        return buf;
    }
    
    if (n == DOWN) {
        if ((my_pos == cLU) || (my_pos == cRU)) {
            return buf;
        }
        TreeIndex c = get_parent();
        buf.push_back(c.get_face_neighbor(n));
        return buf;
    }

    return buf;
}

vector<TreeIndex> TreeIndex::get_larger_possible_corner_neighbour(CornerNeigh n) {
    vector<TreeIndex> ret;

    Child my_pos = get_child_pos();
    if (n == LD) {
        if (my_pos != cRU) {
            return ret;
        }
    }
    if (n == LU) {
        if (my_pos != cRD) {
            return ret;
        }
    }
    if (n == LD) {
        if (my_pos != cRU) {
            return ret;
        }
    }
    if (n == LD) {
        if (my_pos != cRU) {
            return ret;
        }
    }

    TreeIndex c = get_parent();
    ret.push_back(c.get_corner_neighbor(n));
    return ret;
}

vector<TreeIndex> TreeIndex::get_halfsize_possible_face_neighbours(Neigh n) {
    TreeIndex fullSizeNeigh = get_face_neighbor(n);
    vector<TreeIndex> halfSizeNeighs;
    halfSizeNeighs.reserve( 2 );
    if (n == DOWN) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(cLU));
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(cRU));
        return halfSizeNeighs;
    }
    if (n == UP) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(cLD));
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(cRD));
        return halfSizeNeighs;
    }
    if (n == LEFT) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(cRU));
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(cRD));
        return halfSizeNeighs;
    }
    if (n == RIGHT) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(cLU));
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(cLD));
        return halfSizeNeighs;
    }
    return halfSizeNeighs;
}

vector<TreeIndex>& TreeIndex::get_halfsize_possible_face_neighbours_optimized(vector<TreeIndex>& buf, Neigh n) {
    TreeIndex fullSizeNeigh = get_face_neighbor(n);
    if (n == DOWN) {
        buf.push_back(fullSizeNeigh.get_child(cLU));
        buf.push_back(fullSizeNeigh.get_child(cRU));
        return buf;
    }
    if (n == UP) {
        buf.push_back(fullSizeNeigh.get_child(cLD));
        buf.push_back(fullSizeNeigh.get_child(cRD));
        return buf;
    }
    if (n == LEFT) {
        buf.push_back(fullSizeNeigh.get_child(cRU));
        buf.push_back(fullSizeNeigh.get_child(cRD));
        return buf;
    }
    if (n == RIGHT) {
        buf.push_back(fullSizeNeigh.get_child(cLU));
        buf.push_back(fullSizeNeigh.get_child(cLD));
        return buf;
    }

    return buf;
}

vector<TreeIndex> TreeIndex::get_halfsize_possible_corner_neighbours(CornerNeigh n) {
    TreeIndex fullSizeNeigh = get_corner_neighbor(n);
    vector<TreeIndex> halfSizeNeighs;
    if (n == LU) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(cRD));
        return halfSizeNeighs;
    }
    if (n == RU) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(cLD));
        return halfSizeNeighs;
    }
    if (n == LD) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(cRU));
        return halfSizeNeighs;
    }
    if (n == RD) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(cLU));
        return halfSizeNeighs;
    }
    return halfSizeNeighs;
}

vector<TreeIndex> TreeIndex::get_all_halfsize_possible_neighs() {
    vector<TreeIndex> res;
    res.reserve( 8 );

    if (!is_left_border()) {
        res = get_halfsize_possible_face_neighbours_optimized(res, LEFT);
    }
    if (!is_right_border()) {
        res = get_halfsize_possible_face_neighbours_optimized(res, RIGHT);
    }
    if (!is_upper_border()) {
        res = get_halfsize_possible_face_neighbours_optimized(res, UP);
    }
    if (!is_down_border()) {
        res = get_halfsize_possible_face_neighbours_optimized(res, DOWN);
    }

    return res;
}

vector<TreeIndex>& TreeIndex::get_all_halfsize_possible_neighs_optimized(vector<TreeIndex>& buf) {

    if (!is_left_border()) {
        buf = get_halfsize_possible_face_neighbours_optimized(buf, LEFT);
    }
    if (!is_right_border()) {
        buf = get_halfsize_possible_face_neighbours_optimized(buf, RIGHT);
    }
    if (!is_upper_border()) {
        buf = get_halfsize_possible_face_neighbours_optimized(buf, UP);
    }
    if (!is_down_border()) {
        buf = get_halfsize_possible_face_neighbours_optimized(buf, DOWN);
    }

    return buf;
}

vector<TreeIndex> TreeIndex::get_all_samesize_possible_neighs() {
    vector<TreeIndex> res;
    res.reserve( 4 );

    if (!is_left_border()) {
        res.push_back(get_face_neighbor(LEFT));
    }
    if (!is_right_border()) {
        res.push_back(get_face_neighbor(RIGHT));
    }
    if (!is_down_border()) {
        res.push_back(get_face_neighbor(DOWN));
    }
    if (!is_upper_border()) {
        res.push_back(get_face_neighbor(UP));
    }

    return res;
}

vector<TreeIndex>& TreeIndex::get_all_samesize_possible_neighs_optimized(vector<TreeIndex>& buf) {

    if (!is_left_border()) {
        buf.push_back(get_face_neighbor(LEFT));
    }
    if (!is_right_border()) {
        buf.push_back(get_face_neighbor(RIGHT));
    }
    if (!is_down_border()) {
        buf.push_back(get_face_neighbor(DOWN));
    }
    if (!is_upper_border()) {
        buf.push_back(get_face_neighbor(UP));
    }

    return buf;
}

vector<TreeIndex> TreeIndex::get_all_larger_possible_neighs() {
    vector<TreeIndex> res;
    res.reserve(4);

    if (!is_left_border()) {
        res = get_larger_possible_face_neighbour_optimized(res, LEFT);
    }
    if (!is_right_border()) {
        res = get_larger_possible_face_neighbour_optimized(res, RIGHT);
    }
    if (!is_upper_border()) {
        res = get_larger_possible_face_neighbour_optimized(res, UP);
    }
    if (!is_down_border()) {
        res = get_larger_possible_face_neighbour_optimized(res, DOWN);
    }

    return res;
}

vector<TreeIndex>& TreeIndex::get_all_larger_possible_neighs_optimized(vector<TreeIndex>& buf) {

    if (!is_left_border()) {
        buf = get_larger_possible_face_neighbour_optimized(buf, LEFT);
    }
    if (!is_right_border()) {
        buf = get_larger_possible_face_neighbour_optimized(buf, RIGHT);
    }
    if (!is_upper_border()) {
        buf = get_larger_possible_face_neighbour_optimized(buf, UP);
    }
    if (!is_down_border()) {
        buf = get_larger_possible_face_neighbour_optimized(buf, DOWN);
    }

    return buf;
}

// TODO optimize copying
vector<GlobalNumber_t> TreeIndex::get_all_possible_neighbours_ids() {
    vector<TreeIndex> all_neighs;
    all_neighs.reserve(20);

    if (lvl < max_lvl) {
        all_neighs = get_all_halfsize_possible_neighs_optimized(all_neighs);
    }
    all_neighs = get_all_samesize_possible_neighs_optimized(all_neighs);
    if (lvl > 1) {
        all_neighs = get_all_larger_possible_neighs_optimized(all_neighs);
    }

    vector<GlobalNumber_t> all_neighs_ids;
    all_neighs_ids.reserve(all_neighs.size());

    // cout << "get_all_possible_neighbours_ids(" << lvl << "," << i << "," << j << ")= { ";
    for (int jjj = 0; jjj < all_neighs.size(); jjj++) {
        TreeIndex c = all_neighs[jjj];
        // cout << "(" << c.lvl << "," << c.i << "," << c.j << "), ";
        all_neighs_ids.push_back(c.get_global_number());
    }
    // cout << " }" << endl;

    sort(all_neighs_ids.begin(), all_neighs_ids.end());
    vector<GlobalNumber_t>::iterator last = std::unique(all_neighs_ids.begin(), all_neighs_ids.end());
    all_neighs_ids.erase(last, all_neighs_ids.end());

    return all_neighs_ids;
}

bool TreeIndex::is_left_border() {
    TreeIndex c = get_face_neighbor(LEFT);
    return (c.j < 0);
}

bool TreeIndex::is_right_border() {
    TreeIndex c = get_face_neighbor(RIGHT);
    return ((c.j & (-1 << max_lvl)) != 0);
}

bool TreeIndex::is_upper_border() {
    TreeIndex c = get_face_neighbor(UP);
    return ((c.i & (-1 << max_lvl)) != 0);
}

bool TreeIndex::is_down_border() {
    TreeIndex c = get_face_neighbor(DOWN);
    return (c.i < 0);
}

bool TreeIndex::is_left_upper_corner() {
    return is_left_border() && is_upper_border();
}

bool TreeIndex::is_right_upper_corner() {
    return is_right_border() && is_upper_border();
}

bool TreeIndex::is_left_down_corner() {
    return is_left_border() && is_down_border();
}

bool TreeIndex::is_right_down_corner() {
    return is_right_border() && is_down_border();
}

bool TreeIndex::is_border() {
    return is_left_border() || is_right_border() || is_upper_border() || is_down_border();
            // is_left_upper_corner() || is_right_upper_corner() || is_left_down_corner() || is_right_down_corner();
}

// i -- Oy j -- Ox
void TreeIndex::get_corner_coords(double *x, double *y) {
//    double lvl_dx = min_dx * pow(2, max_lvl - lvl);
    *x = min_dx * j;
    *y = min_dx * i;
}


/* * * * * * * * РАБОТА С ЯЧЕЙКАМИ И ДЕРЕВОМ * * * * * * * */

double dist(double x1, double y1, double x2, double y2) {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

// i -- Oy j -- Ox
void Cell::get_spacial_coords(double *x, double *y) {
    double lvl_dx = min_dx * pow(2, max_lvl - lvl);
    *x = min_dx * j + lvl_dx/2; 
    *y = min_dx * i + lvl_dx/2; 
}

void Cell::get_border_cond(char *cond_type, double (**cond_func)(double, double, double)) {
    if (is_left_border()) {
        Area::get_border_cond(Area::LEFT, cond_type, cond_func);
        return;
    }

    if (is_right_border()) {
        Area::get_border_cond(Area::RIGHT, cond_type, cond_func);
        return;
    }

    if (is_down_border()) {
        Area::get_border_cond(Area::DOWN, cond_type, cond_func);
        return;
    }

    if (is_upper_border()) {
        Area::get_border_cond(Area::UP, cond_type, cond_func);
        return;
    }

    *cond_type = -1;
}


double get_lvl_dx(int lvl) {
    return  min_dx * pow(2, max_lvl - lvl);
}



/* * * * * * * * РАБОТА БЛОЧНЫМ ДЕРЕВОМ * * * * * * * */


BlockOfCells::BlockOfCells(int _cells_lvl, int _blk_lvl, int _sz, GlobalNumber_t _i):
        cells_lvl(_cells_lvl), i(_i), idx(_blk_lvl, _i), sz(_sz) {
            // cout << "BlockOfCells(" << _cells_lvl << "," << _blk_lvl << "," << _sz << "," << _i << ")\n";
            refine_marks[0] = 0; refine_marks[1] = 0; refine_marks[2] = 0; refine_marks[3] = 0;
            refine_mark=0;
        }

void BlockOfCells::CreateCells(double (*Temp_func)(double, double)) {

    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {

            double x, y;
            get_spacial_coords(i, j, &x, &y);
            // cout << "get_spacial_coords(blk_i=" << this->i << ",cells_lvl=" << cells_lvl << "," << i << "," << j << ") = (" << x << "," << y << ")\n";  

            SimpleCell c(Temp_func(x, y));
            cells.push_back(c);
        }
    }

    // cout << "CC new_size=" << cells.size() << endl;
}

// i -- Oy j -- Ox
void BlockOfCells::get_spacial_coords(int i, int j, double *x, double *y) {
    double xx, yy;
    idx.get_corner_coords(&xx, &yy);
    // cout << "get_corner_coords(" << this->i << ") = (" << xx << "," << yy << ")\n"; 

    double lvl_dx = min_dx * pow(2, max_lvl - cells_lvl);
    xx += j * lvl_dx + lvl_dx/2;
    yy += i *lvl_dx + lvl_dx/2;

    *x = xx;
    *y = yy;
}


int BlockOfCells::MarkToRefine() {

    if (cells_lvl == max_lvl) {
        return 0;
    }

    int n_of_marks = 0;

    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {

            double x, y;
            get_spacial_coords(i, j, &x, &y);

            if ((cells_lvl == base_lvl) && (Area::Refine1(x, y))) {
                cells[i*sz + j].refine_mark = 1;
//                cells[i*sz + j].temp[0] = 100;
                // cout << int(cells[i].refine_mark)  << " ";
                n_of_marks++;
                mark_quarter(i, j);
            }

            if ((cells_lvl == base_lvl + 1) && (Area::Refine2(x, y))) {
                cells[i*sz + j].refine_mark = 1;
                n_of_marks++;
                mark_quarter(i, j);
            }

            if ((cells_lvl == base_lvl + 2) && (Area::Refine3(x, y))) {
                cells[i*sz + j].refine_mark = 1;
                n_of_marks++;
                mark_quarter(i, j);
            }
        }
    }

    if (n_of_marks > 0) {
        refine_mark = 1;
        return 1;
    }

    return 0;
}

void BlockOfCells::mark_quarter(int i, int j) {
    if ( (i < sz/2) && (j < sz/2) ) {
        refine_marks[0] = 1;
        return;
    }
    if ( (i < sz/2) && (j >= sz/2) ) {
        refine_marks[1] = 1;
        return;
    }
    if ( (i >= sz/2) && (j < sz/2) ) {
        refine_marks[2] = 1;
        return;
    }
    if ( (i >= sz/2) && (j >= sz/2) ) {
        refine_marks[3] = 1;
        return;
    }
}

vector<BlockOfCells> BlockOfCells::Split(double (*Temp_func)(double, double)) {
    vector<BlockOfCells> children;

//? здесь бага ( sz/2 = 0 ) если он был равен 1 (т.е. уровень блока был равен уровню ячейки)

    BlockOfCells bc00 = BlockOfCells(cells_lvl, idx.lvl+1, sz/2, idx.get_child(cLD).get_global_number());
    bc00.CreateCells(Temp_func);
    if (refine_marks[0]) {
        bc00.refine_mark = 1;
    }
    children.push_back(bc00);

    BlockOfCells bc01 = BlockOfCells(cells_lvl, idx.lvl+1, sz/2, idx.get_child(cRD).get_global_number());
    bc01.CreateCells(Temp_func);
    if (refine_marks[1]) {
        bc01.refine_mark = 1;
    }
    children.push_back(bc01);

    BlockOfCells bc10 = BlockOfCells(cells_lvl, idx.lvl+1, sz/2, idx.get_child(cLU).get_global_number());
    bc10.CreateCells(Temp_func);
    if (refine_marks[2]) {
        bc10.refine_mark = 1;
    }
    children.push_back(bc10);

    BlockOfCells bc11 = BlockOfCells(cells_lvl, idx.lvl+1, sz/2, idx.get_child(cRU).get_global_number());
    bc11.CreateCells(Temp_func);
    if (refine_marks[3]) {
        bc11.refine_mark = 1;
    }
    children.push_back(bc11);

    return children;
}

void BlockOfCells::RefineCells() {
    cells = vector<SimpleCell>();

    cells_lvl++;
    sz *= 2;
    CreateCells(&Area::T0);
}

void BlockOfCells::ClearMarks() {
    refine_mark = 0;
    refine_marks[0] = 0; refine_marks[1] = 0; refine_marks[2] = 0; refine_marks[3] = 0;
}



BlockedLinearTree::BlockedLinearTree(double (*Temp_func)(double, double)) {
    cout << "building base blocked tree\n";
    max_present_lvl = base_lvl;
    max_present_blk_lvl = base_blk_lvl;

    int global_i_start = 0;
    int global_i_stop  = 1 << (max_lvl*2);
    int global_i_step  = 1 << (2*max_lvl - 2*base_blk_lvl);


    for (int i = global_i_start; i < global_i_stop; i += global_i_step) {
        BlockOfCells bc = BlockOfCells(base_lvl, base_blk_lvl, base_blk_sz, i);
        bc.CreateCells(Temp_func);
        blocks.push_back(bc);
    }
    cout << "base blocked tree built, n_of_blocks=" << blocks.size() << endl;
}

int BlockedLinearTree::MarkToRefine() {
    cout << "analysing refine\n";
    if (max_present_lvl == max_lvl) {
        return 0;
    }

    int n_of_block_marks = 0;
    for (int i = 0; i < blocks.size(); i++) {
        n_of_block_marks += blocks[i].MarkToRefine();
    }

    // cout << "refine analysed, n_of_block_marks=" << n_of_block_marks << " block_marks={";
    // for (BlockOfCells& blk: blocks) {
    //     if (blk.refine_mark) {
    //         cout << blk.idx.get_global_number() << " ";
    //     }
    // }
    // cout << "}\n";
    return (n_of_block_marks > 0);
}

void BlockedLinearTree::RefineBlocks() {
    cout << "refining blocks\n";

    if (max_present_blk_lvl == max_blk_lvl) {
        cout << "max blk level reached, not refining blocks\n";
        return;
    }

    int new_max_present_blk_lvl = max_present_blk_lvl;
    for (int i = 0; i < blocks.size();) {
        // cout << i << " ";
        if (blocks[i].refine_mark) {
            if (blocks[i].idx.lvl < max_blk_lvl) {
                if (blocks[i].idx.lvl + 1 > new_max_present_blk_lvl) {
                    new_max_present_blk_lvl = blocks[i].idx.lvl + 1;
                }
                vector<BlockOfCells> new_blocks = blocks[i].Split(&Area::T0);
                blocks.erase(blocks.begin()+i);
                blocks.insert(blocks.begin()+(i), new_blocks.begin(), new_blocks.end());
                i+= 4;

            } else {
                i++;
            }
        } else {
            i++;
        }
    }

    max_present_blk_lvl = new_max_present_blk_lvl;
    cout << "blocks refined, new max_present_blk_lvl = " << max_present_blk_lvl << endl;
}

void BlockedLinearTree::RefineCells() {
    cout << "refining cells\n";
    int new_max_present_lvl = max_present_lvl;

    for (int i = 0; i < blocks.size(); i++) {
        if (blocks[i].refine_mark) {
            // cout << "refinigh cells in blk_idx=" << blocks[i].idx.get_global_number() << endl;
            blocks[i].RefineCells();
            if (blocks[i].cells_lvl > new_max_present_lvl) {
                new_max_present_lvl = blocks[i].cells_lvl;
            }
            blocks[i].ClearMarks();
        }
    }

    max_present_lvl = new_max_present_lvl;
    cout << "cells refined, new max_present_lvl = " << max_present_lvl << endl;
}

void BlockedLinearTree::BuildNeighs() {
    cout << "building neighs\n";

    for (int i = 0; i < blocks.size(); i++) {

        // cout << "BN " << i << " h1\n";
        if (!blocks[i].idx.is_left_border()) {
            // cout << "BN " << i << " h11\n";
            TreeIndex possible_neigh_idx = blocks[i].idx.get_face_neighbor(LEFT);
            // cout << "BN " << i << " h12\n";
            BlockOfCells* possible_neigh = find_block(possible_neigh_idx);
            // cout << "BN " << i << " h13 idx=" << possible_neigh_idx.get_global_number() << std::endl;
            if (possible_neigh != NULL) {
                blocks[i].neighs_left = find_block_children(possible_neigh_idx, RIGHT);
            }
            // cout << "BN " << i << " h14\n";
        }

        // cout << "BN " << i << " h2\n";
        if (!blocks[i].idx.is_right_border()) {
            TreeIndex possible_neigh_idx = blocks[i].idx.get_face_neighbor(RIGHT);
            BlockOfCells* possible_neigh = find_block(possible_neigh_idx);
            if (possible_neigh != NULL) {
                blocks[i].neighs_right = find_block_children(possible_neigh_idx, LEFT);
            }
        }

        // cout << "BN " << i << " h3\n";
        if (!blocks[i].idx.is_upper_border()) {
            TreeIndex possible_neigh_idx = blocks[i].idx.get_face_neighbor(UP);
            BlockOfCells* possible_neigh = find_block(possible_neigh_idx);
            if (possible_neigh != NULL) {
                blocks[i].neighs_upper = find_block_children(possible_neigh_idx, DOWN);
            }
        }

        // cout << "BN " << i << " h3\n";
        if (!blocks[i].idx.is_down_border()) {
            TreeIndex possible_neigh_idx = blocks[i].idx.get_face_neighbor(DOWN);
            // if (blocks[i].idx.get_global_number() == GlobalNumber_t(2048)) {
            //     cout << "For 2048 down: " << possible_neigh_idx.get_global_number() << " " << possible_neigh_idx.lvl << endl;
            // }
            BlockOfCells* possible_neigh = find_block(possible_neigh_idx);
            // if (blocks[i].idx.get_global_number() == GlobalNumber_t(2048)) {
            //     cout << "find(512) = " << possible_neigh << endl;
            // }
            if (possible_neigh != NULL) {
                blocks[i].neighs_down = find_block_children(possible_neigh_idx, UP);
            }
        }

    }

    // for (BlockOfCells& blk: blocks) {
    //     // cout << "h1\n";
    //     cout << blk.idx.get_global_number() << " {";
    //     for (BlockOfCells *b_ptr: blk.neighs_down) {
    //         cout << b_ptr->idx.get_global_number() << ",";
    //     }
    //     cout << "} ";
    //     // cout << "h1\n";
    //     cout << blk.idx.get_global_number() << " {";
    //     for (BlockOfCells *b_ptr: blk.neighs_upper) {
    //         cout << b_ptr->idx.get_global_number() << ",";
    //     }
    //     cout << "} ";
    //     cout << blk.idx.get_global_number() << " {";
    //     for (BlockOfCells *b_ptr: blk.neighs_left) {
    //         cout << b_ptr->idx.get_global_number() << ",";
    //     }
    //     cout << "} ";
    //     cout << blk.idx.get_global_number() << " {";
    //     for (BlockOfCells *b_ptr: blk.neighs_right) {
    //         cout << b_ptr->idx.get_global_number() << ",";
    //     }
    //     cout << "} \n";
    // }

    // cout << "BN " << " h4\n";
    for (int i = 0; i < blocks.size(); i++) {
        // left
        for (int jjj = 0; jjj < blocks[i].neighs_left.size(); jjj++) {
            BlockOfCells *b_ptr = blocks[i].neighs_left[jjj];
            if(std::find(b_ptr->neighs_right.begin(), b_ptr->neighs_right.end(), &blocks[i]) == b_ptr->neighs_right.end()) {
                // не делаю сортировку потому что в таком случае должен
                // быть всего один сосед больший по размеру
                b_ptr->neighs_right.insert(b_ptr->neighs_right.begin(), &blocks[i]);
            }
        }

        // right
        for (int jjj = 0; jjj < blocks[i].neighs_right.size(); jjj++) {
            BlockOfCells *b_ptr = blocks[i].neighs_right[jjj];
            if(std::find(b_ptr->neighs_left.begin(), b_ptr->neighs_left.end(), &blocks[i]) == b_ptr->neighs_left.end()) {
                b_ptr->neighs_left.insert(b_ptr->neighs_left.begin(), &blocks[i]);
            }
        }

        // upper
        for (int jjj = 0; jjj < blocks[i].neighs_upper.size(); jjj++) {
            BlockOfCells *b_ptr = blocks[i].neighs_upper[jjj];
            if(std::find(b_ptr->neighs_down.begin(), b_ptr->neighs_down.end(), &blocks[i]) == b_ptr->neighs_down.end()) {
                b_ptr->neighs_down.insert(b_ptr->neighs_down.begin(), &blocks[i]);
            }
        }

        // down
        for (int jjj = 0; jjj < blocks[i].neighs_down.size(); jjj++) {
            BlockOfCells *b_ptr = blocks[i].neighs_down[jjj];
            if(std::find(b_ptr->neighs_upper.begin(), b_ptr->neighs_upper.end(), &blocks[i]) == b_ptr->neighs_upper.end()) {
                b_ptr->neighs_upper.insert(b_ptr->neighs_upper.begin(), &blocks[i]);
            }
        }
    }

    // for (BlockOfCells& blk: blocks) {
    //     cout << blk.idx.get_global_number() << " {";
    //     for (BlockOfCells *b_ptr: blk.neighs_down) {
    //         cout << b_ptr->idx.get_global_number() << ",";
    //     }
    //     cout << "} ";
    //     cout << blk.idx.get_global_number() << " {";
    //     for (BlockOfCells *b_ptr: blk.neighs_upper) {
    //         cout << b_ptr->idx.get_global_number() << ",";
    //     }
    //     cout << "} ";
    //     cout << blk.idx.get_global_number() << " {";
    //     for (BlockOfCells *b_ptr: blk.neighs_left) {
    //         cout << b_ptr->idx.get_global_number() << ",";
    //     }
    //     cout << "} ";
    //     cout << blk.idx.get_global_number() << " {";
    //     for (BlockOfCells *b_ptr: blk.neighs_right) {
    //         cout << b_ptr->idx.get_global_number() << ",";
    //     }
    //     cout << "} \n";
    // }

    cout << "built neighs\n";
}

BlockOfCells* BlockedLinearTree::find_block(TreeIndex target) {

    GlobalNumber_t t = target.get_global_number();

    // binary search
    int left = 0;
    int right = blocks.size();
    while (left < right) {
        int midi = (right + left) / 2;
        BlockOfCells& mid = blocks[midi];
        GlobalNumber_t val = mid.idx.get_global_number();
        if (val == t) {
            return &mid;
        }
        if (val < t) {
            left = midi + 1;
        } else {
            right = midi;
        }
    }
    return NULL;
}

BlockOfCells* BlockedLinearTree::find_block(GlobalNumber_t target) {

    GlobalNumber_t t = target;

    // binary search
    int left = 0;
    int right = blocks.size();
    while (left < right) {
        int midi = (right + left) / 2;
        BlockOfCells& mid = blocks[midi];
        GlobalNumber_t val = mid.idx.get_global_number();
        if (val == t) {
            return &mid;
        }
        if (val < t) {
            left = midi + 1;
        } else {
            right = midi;
        }
    }
    return NULL;
}

vector<BlockOfCells*> BlockedLinearTree::find_block_children(TreeIndex target, Neigh n) {

    // cout << "find_block_children( (" << target.lvl << "," << target.i << "," << target.j << "), " << n <<  ")" << endl;
    // cout << "find_block_children(" << target.get_global_number() << ")" << endl;

    vector<TreeIndex> res_idxs;
    res_idxs.push_back(target);

    int changed = 1;
    while (changed) {
        changed = 0;

        int i = 0;
        while (i < res_idxs.size()) {
            TreeIndex cur_idx = res_idxs[i];

            // cout << res_list[i].idx << endl;
            if (cur_idx.lvl == max_present_blk_lvl) {
                i++;
                continue;
            }

            // if (cur_idx.get_global_number() == GlobalNumber_t(3328)) {
            //     for (TreeIndex ii: res_idxs) {
            //         cout << "  " << ii.get_global_number();
            //     }
            //     cout << endl;
            // }

            // find two needed children
            Child ch1, ch2;
            if (n == LEFT) {
                ch1 = cLD; ch2 = cLU;
            } else if (n == RIGHT) {
                ch1 = cRD; ch2 = cRU;
            } else if (n == UP) {
                ch1 = cLU; ch2 = cRU;
            } else if (n == DOWN) {
                ch1 = cLD; ch2 = cRD;
            }

            GlobalNumber_t t1 = cur_idx.get_child(ch1).get_global_number();
            GlobalNumber_t t2 = cur_idx.get_child(ch2).get_global_number();

            BlockOfCells *b1 = NULL, *b2 = NULL;
            b1 = find_block(t1);
            b2 = find_block(t2);

            // if (target.get_global_number() == GlobalNumber_t(512)) {
            //     cout << "For 512: " << cur_idx.get_global_number() << " " << t1 << " " << t2 << ";" << b1 << " " << b2 << endl;
            // }

            // update children list
            if (b1 && b2) {
                res_idxs.erase(res_idxs.begin()+i);
                res_idxs.insert(res_idxs.begin()+i, cur_idx.get_child(ch2));
                res_idxs.insert(res_idxs.begin()+i, cur_idx.get_child(ch1));

                changed = 1;
                i += 2;
            } else {
                i++;
            }
        }
    }

    vector<BlockOfCells*> res;
    for (int jjj = 0; jjj < res_idxs.size(); jjj++) {
        // cout << "1";
        res.push_back(find_block(res_idxs[jjj]));
    }
    // cout << res.size() << endl;

    return res;
}

void BlockedLinearTree::Decompose(int n_procs) {
    cout << "decomposing grid\n";

    int sum_weight = 0;
    for (int i = 0; i < blocks.size(); i++) {
        sum_weight += blocks[i].GetNOfCells();
    }

    cout << "N_OF_BLOCKS=" << blocks.size() << endl;
    cout << "N_OF_CELLS=" << sum_weight << endl;

    proc_blocks = vector<int>(blocks.size(), 0);

    int optimal_weight = sum_weight / n_procs;
    cout << "OPTIMAL_WEIGHT=" << optimal_weight << endl;
    
    int max_weight = 0;
    int max_weight_proc = -1;
    int cur_weight = 0;
    int proc_i = 0;

    vector<int> weights(n_procs, 0);

    for (int i = 0; i < blocks.size(); i++) {
        proc_blocks[i] = proc_i;
        cur_weight += blocks[i].GetNOfCells();
        weights[proc_i] += blocks[i].GetNOfCells();
        if (cur_weight >= optimal_weight) {
            if (cur_weight > max_weight) {
                max_weight = cur_weight;
                max_weight_proc = proc_i;
            }
            proc_i++;
            cur_weight = 0;
        }
    }

    cout << "WEIGHTS = [ ";
    for (int i = 0; i < n_procs; i++) {
        cout << weights[i] << " ";
    }
    cout << "]" << endl;

    // for (int i = 0; i < proc_blocks.size(); i++) {
    //     cout << proc_blocks[i] << " ";
    // }
    // cout << endl;

    cout << "MAX_WEIGHT=" << max_weight << " MAX_WEIGHT_PROC=" << max_weight_proc << endl;

    cout << "grid decomposed\n";
}

void BlockedLinearTree::Write(string filename) {
    cout << "writing grid\n";
    vector<char> buf = GenWriteStruct(0);
    int len = buf.size();

    std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary);
    fout.write(&buf[0], len);
    fout.close();
    cout << "grid written\n";
}

vector<char> BlockedLinearTree::GenWriteStruct(char light) {
    vector<char> buf;
    // std::cout << " gen structs for write start\n";
    for (int blk_i = 0; blk_i < blocks.size(); blk_i++) {
        int lvl = blocks[blk_i].sz / 2;
        for (int k = 0; k < blocks[blk_i].GetNOfCells(); k++) {
            int cell_i = 0, cell_j = 0;
            split_ints(k, lvl, &cell_i, &cell_j);
            SimpleCell c = blocks[blk_i].cells[cell_i * blocks[blk_i].sz + cell_j];

            int _lvl = blocks[blk_i].cells_lvl;
            char *tmp = (char *)(&_lvl);
            for (int i = 0; i < sizeof(int); i++) {
                buf.push_back(tmp[i]);
            }

            int tree_i, tree_j;
            GlobalNumber_t glob_idx = get_glob_idx(blocks[blk_i].i, k, _lvl);
            split_ints(glob_idx, &tree_i, &tree_j);
            // cout << "blk_i=" << blocks[blk_i].i << " cell_k=" << k << " cell_i=" << cell_i << " cell_j=" << cell_j << " tree_i=" << tree_i << " tree_j=" << tree_j << endl;

            tmp = (char *)(&tree_i);
            for (int i = 0; i < sizeof(int); i++) {
                buf.push_back(tmp[i]);
            }
            tmp = (char *)(&tree_j);
            for (int i = 0; i < sizeof(int); i++) {
                buf.push_back(tmp[i]);
            }

            int blk_n = blk_i;
            tmp = (char *)(&blk_n);
            for (int i = 0; i < sizeof(int); i++) {
                buf.push_back(tmp[i]);
            }

            int proc_n = 0;
            if (!light) {
                proc_n = proc_blocks[blk_i];
            }
            tmp = (char *)(&proc_n);
            for (int i = 0; i < sizeof(int); i++) {
                buf.push_back(tmp[i]);
            }

            tmp = (char *)(&c.temp[0]);
            for (int i = 0; i < sizeof(double); i++) {
                buf.push_back(tmp[i]);
            }
        }
    }
    // std::cout << " gen structs for write finished\n";
    return buf;
}

SimpleCell *BlockOfCells::find_border_cell_by_global_idx(GlobalNumber_t target) {
    SimpleCell *res = NULL;

    GlobalNumber_t mask(-1);
    mask = mask << 2*(max_lvl - idx.lvl); // тут наверное аналогично check_cell_owner сдвиг
    mask = ~mask;

    GlobalNumber_t idx = target & mask;
    idx = idx >> (2*(max_lvl - cells_lvl));
    if ((idx > sz*sz) || (idx < 0)) {
        cout << "[ERROR] idx out of bounds\n";
        return NULL;
    }

    // по идее тут нужен сплит интс
    int i, j;
    split_ints(idx, max_lvl, &i, &j);

    return &cells[i*sz+j]; // тут нужен индекс в массиве ячеек (который тип матрица)
}

// get_glob_idx получает глобальный номер ячейки во всей сетке по глобальному 
// номеру блока и глобавльному номеру ячейки в блоке (+ уровень ячеек)
GlobalNumber_t get_glob_idx(GlobalNumber_t blk_i, GlobalNumber_t cell_i, int cell_lvl) {
    GlobalNumber_t res = blk_i;
    // res = res << 2*(max_lvl - max_blk_lvl);
    res = res | (cell_i << 2*(max_lvl - cell_lvl));
    return res;
}

GlobalNumber_t get_cell_global_number_in_full_grid(GlobalNumber_t blk_gloal_number, int cell_i, int cell_j, int cell_lvl) {
    TreeIndex cell_idx = TreeIndex(cell_lvl, cell_i, cell_j);
    GlobalNumber_t cell_global_number = cell_idx.get_global_number();
    return get_glob_idx(blk_gloal_number, cell_global_number, cell_lvl);
}

int check_cell_owner(GlobalNumber_t blk_i, int blk_lvl, GlobalNumber_t cell_i) {
    GlobalNumber_t mask(-1);

    // std::bitset<sizeof(GlobalNumber_t)*8> x(mask);
    // cout << "sizeof(GlobalNumber_t)=" << sizeof(GlobalNumber_t) << " sizeof(long long int)=" << sizeof(long long int) << " mask:" << x << '\n';

    mask = mask << 2*(max_lvl - blk_lvl);

    // std::bitset<sizeof(GlobalNumber_t)*8> bs_cell_i(cell_i);
    // std::bitset<sizeof(GlobalNumber_t)*8> bs_blk_i(blk_i);
    // std::bitset<sizeof(GlobalNumber_t)*8> bs_mask(mask);
    // GlobalNumber_t cell_blk_i = cell_i & mask;
    // std::bitset<sizeof(GlobalNumber_t)*8> bs_cell_blk_i(cell_blk_i);
    // cout << "check_cell_owner: " << bs_cell_i << " & " << bs_mask << " = " << bs_cell_blk_i << " ?== " << bs_blk_i << '\n';

    return (cell_i & mask) ==  blk_i;
}

// НЕ ОБРАЩАТЬСЯ К ЯЧЕЙКАМ БЛОКА В ЭТОЙ ФУНКЦИИ!!!
vector<GlobalNumber_t> find_cell_neighs_ids_in_blk(GlobalNumber_t cell_glob_idx, int cell_lvl, BlockOfCells* blk, Neigh neigh_dir) {

    int neigh_lvl = blk->cells_lvl;

    // if (neigh_dir == DOWN) {
    //     cout << "find_cell_neighs_ids_in_blk DOWN: cell_idx=" << cell_glob_idx << " cell_lvl=" << cell_lvl << " neigh_blk_idx=" << blk->idx.get_global_number() << " neigh_lvl=" << blk->cells_lvl << endl;
    // }

    if (abs(cell_lvl - neigh_lvl) > 1) {
        cout << "[ERROR] not 2:1 balanced grid!\n";
        exit(-1);
    }

    TreeIndex c_idx = TreeIndex(cell_lvl, cell_glob_idx);
    vector<TreeIndex> neighs;
    if (neigh_lvl > cell_lvl) {
        neighs = c_idx.get_halfsize_possible_face_neighbours(neigh_dir);
    } else if (neigh_lvl < cell_lvl) {
        neighs = c_idx.get_larger_possible_face_neighbour(neigh_dir);
    } else {
        neighs.push_back(c_idx.get_face_neighbor(neigh_dir));
    }

    // if (neigh_dir == DOWN) {
    //     cout << "possible neighs: ";
    //     for (auto nei: neighs) {
    //         cout << nei.get_global_number() << " ";
    //     }
    //     cout << endl;
    // }

    vector<GlobalNumber_t> res;
    for (int jjj = 0; jjj < neighs.size(); jjj++) {
        GlobalNumber_t n = neighs[jjj].get_global_number();
        if (check_cell_owner(blk->idx.get_global_number(), blk->idx.lvl, n)) {
            res.push_back(n);
        }
    }

    // if (neigh_dir == DOWN) {
    //     cout << "real neighs: ";
    //     for (auto nei: res) {
    //         cout << nei << " ";
    //     }
    //     cout << endl;
    // }

    return res;
}


void BlockedLinearTree::WriteBlocks(string filename) {
    vector<char> buf = GenWriteBlocksStruct();

    int len = buf.size();
    cout << "WriteBlocks len=" << len << endl;

    std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary);
    fout.write(&buf[0], len);
    fout.close();
}



vector<char> BlockedLinearTree::GenWriteBlocksStruct() {
    vector<char> buf;
    std::cout << " gen structs for write start\n";

    for (int blk_i = 0; blk_i < blocks.size(); blk_i++) {
        
        int _lvl = blocks[blk_i].idx.lvl;
        char *tmp = (char *)(&_lvl);
        for (int i = 0; i < sizeof(int); i++) {
            buf.push_back(tmp[i]);
        }

        int _i = blocks[blk_i].idx.i;
        tmp = (char *)(&_i);
        for (int i = 0; i < sizeof(int); i++) {
            buf.push_back(tmp[i]);
        }

        int _j = blocks[blk_i].idx.j;
        tmp = (char *)(&_j);
        for (int i = 0; i < sizeof(int); i++) {
            buf.push_back(tmp[i]);
        }

        int _c_lvl = blocks[blk_i].cells_lvl;
        tmp = (char *)(&_c_lvl);
        for (int i = 0; i < sizeof(int); i++) {
            buf.push_back(tmp[i]);
        }

        int _sz = blocks[blk_i].sz;
        tmp = (char *)(&_sz);
        for (int i = 0; i < sizeof(int); i++) {
            buf.push_back(tmp[i]);
        }

        // cout << "writing blk(" << _lvl << " " << _i << " " << _j << " " << _c_lvl << " " << _sz << ")" << endl;

        // left neighs
        for (int i = 0; i < pow(2, max_blk_lvl-base_blk_lvl); i++) {
            GlobalNumber_t _neigh(-1);
            if (i < blocks[blk_i].neighs_left.size()) {
                _neigh = blocks[blk_i].neighs_left[i]->idx.get_global_number();
            }
            tmp = (char *)(&_neigh);
            for (int j = 0; j < sizeof(GlobalNumber_t); j++) {
                buf.push_back(tmp[j]);
            }
        }

        // right neighs
        for (int i = 0; i < pow(2, max_blk_lvl-base_blk_lvl); i++) {
            GlobalNumber_t _neigh(-1);
            if (i < blocks[blk_i].neighs_right.size()) {
                _neigh = blocks[blk_i].neighs_right[i]->idx.get_global_number();
            }
            tmp = (char *)(&_neigh);
            for (int j = 0; j < sizeof(GlobalNumber_t); j++) {
                buf.push_back(tmp[j]);
            }
        }

        // upper neighs
        for (int i = 0; i < pow(2, max_blk_lvl-base_blk_lvl); i++) {
            GlobalNumber_t _neigh(-1);
            if (i < blocks[blk_i].neighs_upper.size()) {
                _neigh = blocks[blk_i].neighs_upper[i]->idx.get_global_number();
            }
            tmp = (char *)(&_neigh);
            for (int j = 0; j < sizeof(GlobalNumber_t); j++) {
                buf.push_back(tmp[j]);
            }
        }

        // down neighs
        for (int i = 0; i < pow(2, max_blk_lvl-base_blk_lvl); i++) {
            GlobalNumber_t _neigh(-1);
            if (i < blocks[blk_i].neighs_down.size()) {
                _neigh = blocks[blk_i].neighs_down[i]->idx.get_global_number();
            }
            // cout << "***writing down_neigh " << _neigh << endl;
            tmp = (char *)(&_neigh);
            for (int j = 0; j < sizeof(GlobalNumber_t); j++) {
                buf.push_back(tmp[j]);
            }
        }
    }

    std::cout << " gen structs for write finished\n";
    return buf;
}

void BlockedLinearTree::WriteOffsets(string filename) {
    int rec_len = 5 * sizeof(int) + 4 * pow(2, max_blk_lvl-base_blk_lvl) * sizeof(GlobalNumber_t);

    vector<int> offsets_lens;

    int cur_proc = 0;
    int cur_len = 0;
    int cur_offset = 0;
    offsets_lens.push_back(0);
    for (int i = 0; i < proc_blocks.size(); i++) {
        // cout << "proc block " << i  << " proc=" << proc_blocks[i] << endl;
        if (proc_blocks[i] != cur_proc) {
            offsets_lens.push_back(cur_len*rec_len);
            cur_offset += cur_len;
            offsets_lens.push_back(cur_offset*rec_len);
            cur_len = 0;
            cur_proc++;
        }
        cur_len++;
    }
    offsets_lens.push_back(cur_len*rec_len);
    // offsets_lens = vector<int>(offsets_lens.begin()+2, offsets_lens.end());

    cout << "OFFSETS={ ";
    for (int i = 0; i < offsets_lens.size(); i++) {
        cout << offsets_lens[i] << " ";
    }
    cout << "}" << endl;

    std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary);
    fout.write((char *)&offsets_lens[0], offsets_lens.size() * sizeof(int));
    fout.close();
}


void BlockedLinearTree::GenFromWriteBlocksStruct(vector<char> buf, double (*start_func)(double, double)) {

    max_present_lvl = base_lvl;
    max_present_blk_lvl = base_blk_lvl;

    int n_neighs = pow(2, max_blk_lvl-base_blk_lvl);
    cout << "n_neighs=" << n_neighs << endl;
    int one_sz = 5 * sizeof(int) + 4 * n_neighs * sizeof(GlobalNumber_t);

    char *p = &buf[0];
    int pos = 0;

    int lvl_offset   = 0;
    int i_offset     = sizeof(int);
    int j_offset     = 2 * sizeof(int);
    int c_lvl_offset = 3 * sizeof(int);
    int sz_offset    = 4 * sizeof(int);
    int left_offset  = sz_offset    + sizeof(int);
    int right_offset = left_offset  + n_neighs * sizeof(GlobalNumber_t);
    int upper_offset = right_offset + n_neighs * sizeof(GlobalNumber_t);
    int down_offset  = upper_offset + n_neighs * sizeof(GlobalNumber_t);


    while (pos < buf.size()) {
        TreeIndex idx;
        idx.lvl  = * ((int *)(&p[pos+lvl_offset]));
        idx.i    = * ((int *)(&p[pos+i_offset]));
        idx.j    = * ((int *)(&p[pos+j_offset]));

        int cells_lvl = * ((int *)(&p[pos+c_lvl_offset]));
        int sz = * ((int *)(&p[pos+sz_offset]));

        BlockOfCells blk(idx, cells_lvl, sz);
        blk.CreateCells(start_func);

        cout << blk.idx.get_global_number() << "  ";

        GlobalNumber_t tmp(0);
        // cout << "reading left neighs... ";
        for (int j = 0; j < n_neighs; j++) {
            tmp = * ((GlobalNumber_t *)(&p[pos+left_offset+j* sizeof(GlobalNumber_t)]));
            // cout << tmp << " ";
            if (tmp != -1) {
                blk.neighs_left_idxs.push_back(tmp);
            }
        }
        // cout << "|";

        // cout << "reading right neighs... ";
        for (int j = 0; j < n_neighs; j++) {
            tmp = * ((GlobalNumber_t *)(&p[pos+right_offset+j* sizeof(GlobalNumber_t)]));
            // cout << tmp << " ";
            if (tmp != -1) {
                blk.neighs_right_idxs.push_back(tmp);
            }
        }
        // cout << endl;

        for (int j = 0; j < n_neighs; j++) {
            tmp = * ((GlobalNumber_t *)(&p[pos+upper_offset+j* sizeof(GlobalNumber_t)]));
            if (tmp != -1) {
                blk.neighs_upper_idxs.push_back(tmp);
            }
        }

        // cout << "reading down neighs... ";
        for (int j = 0; j < n_neighs; j++) {
            tmp = * ((GlobalNumber_t *)(&p[pos+down_offset+j* sizeof(GlobalNumber_t)]));
            // cout << tmp << " ";
            if (tmp != -1) {
                blk.neighs_down_idxs.push_back(tmp);
            }
        }
        // cout << endl;

        // cout << " BlockOfCells created {" << blk.idx.get_global_number() << ", " 
        //     << blk.cells_lvl << "," <<  blk.sz << "; " << blk.cells.size() << "; " << blk.neighs_down_idxs.size() << "}\n";

        blocks.push_back(blk);

        // std::cout << "cell pushed pos=" << pos << endl;
        if (idx.lvl > max_present_blk_lvl) {
            max_present_blk_lvl = idx.lvl;
        }

        if (idx.lvl > max_present_blk_lvl) {
            max_present_blk_lvl = idx.lvl;
        }

        pos += one_sz;
    }

    // cout << "first cell = Cell(" << cells[0].lvl << ", " << cells[0].i << "," << cells[0].j << ", " << cells[0].temp[0] << ")";
    // cout << "last cell = Cell(" << cells[cells.size()-1].lvl << ", " << cells[cells.size()-1].i << "," << cells[cells.size()-1].j << ", " << cells[cells.size()-1].temp[0] << ")";
}


//double BlockedLinearTree::BuildBlocks(int block_size) {
//    // анализируем возможность объединения последовательных ячеек
//    // TODO нужна структура сохраняющая уровни блоков кроме начал
//
//    map<int, BlockOfCells> lvl_blocks;
//
//    // base grid
//    if (base_sz%block_size != 0)  {
//        cout << "base_size=" << base_sz << " not divisible into block_sz=" << block_sz << endl;
//    }
//    int base_blocks_n = base_sz / block_size;
//    lvl_blocks[base_lvl] = vector<vector<int>>()
//
//    vector<int> seq_lvl_offsets;
//
//    cur_lvl = 0;
//    cur_block_sz = 0;
//    for (int i = 0; i < cells.size(); i++) {
//        if (cells[i].lvl != cur_lvl) {
//            // mark one level stop
//            seq_lvl_offsets.push_back(i);
//
//            cur_block_sz = 0;
//            cur_lvl = lvl;
//        }
//        cur_block_sz++;
//    }
//
//    for (int lvl_offset: seq_lvl_offsets) {
//
//    }
//
//
//}

