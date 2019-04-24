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
    std::cout << "new TreeIndex(" << lvl << ", " << globalNumber << ")\n";
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

    if (c == Child::cLD) {        // (0,0)
        child.i = i;
        child.j = j;
    } else if (c == Child::cRD) { // (0,1)
        child.i = i;
        child.j = j + h;
    } else if (c == Child::cLU) { // (1,0)
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
        return Child::cLD;
    }
    if ((hi == 0) && (hj != 0)) {
        return Child::cRD;
    }
    if ((hi != 0) && (hj == 0)) {
        return Child::cLU;
    }
    if ((hi != 0) && (hj != 0)) {
        return Child::cRU;
    }
    return Child::cLD;
}

TreeIndex TreeIndex::get_face_neighbor(Neigh n) {
    int h = 1 << (max_lvl - lvl); // 2^(b-l)
    TreeIndex c;
    c.lvl = lvl;
    c.i   = i + ((n == Neigh::DOWN) ? -h : (n == Neigh::UP) ? h : 0);
    c.j   = j + ((n == Neigh::LEFT) ? -h : (n == Neigh::RIGHT) ? h : 0);
    return c;
}

TreeIndex TreeIndex::get_corner_neighbor(CornerNeigh n) {
    int h = 1 << (max_lvl - lvl - 1); // 2^(b-l)
    TreeIndex c;
    c.lvl = lvl;
    c.i   = i + ( ((n == CornerNeigh::LU) || (n == CornerNeigh::LD)) ? -h : h);
    c.j   = j + ( ((n == CornerNeigh::LD) || (n == CornerNeigh::RD)) ? -h : h);
    return c;
}

vector<TreeIndex> TreeIndex::get_larger_possible_face_neighbour(Neigh n) {
    vector<TreeIndex> ret;

    Child my_pos = get_child_pos();
    if (n == Neigh::RIGHT) {
        if ((my_pos == Child::cLD) || (my_pos == Child::cLU)) {
            return ret;
        }
        TreeIndex c = get_parent();
        ret.push_back(c.get_face_neighbor(n));
        return ret;
    }

    if (n == Neigh::LEFT) {
        if ((my_pos == Child::cRD) || (my_pos == Child::cRU)) {
            return ret;
        }
        TreeIndex c = get_parent();
        ret.push_back(c.get_face_neighbor(n));
        return ret;
    }

    if (n == Neigh::UP) {
        if ((my_pos == Child::cLD) || (my_pos == Child::cRD)) {
            return ret;
        }
        TreeIndex c = get_parent();
        ret.push_back(c.get_face_neighbor(n));
        return ret;
    }
    
    if (n == Neigh::DOWN) {
        if ((my_pos == Child::cLU) || (my_pos == Child::cRU)) {
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
    if (n == Neigh::RIGHT) {
        if ((my_pos == Child::cLD) || (my_pos == Child::cLU)) {
            return buf;
        }
        TreeIndex c = get_parent();
        buf.push_back(c.get_face_neighbor(n));
        return buf;
    }

    if (n == Neigh::LEFT) {
        if ((my_pos == Child::cRD) || (my_pos == Child::cRU)) {
            return buf;
        }
        TreeIndex c = get_parent();
        buf.push_back(c.get_face_neighbor(n));
        return buf;
    }

    if (n == Neigh::UP) {
        if ((my_pos == Child::cLD) || (my_pos == Child::cRD)) {
            return buf;
        }
        TreeIndex c = get_parent();
        buf.push_back(c.get_face_neighbor(n));
        return buf;
    }
    
    if (n == Neigh::DOWN) {
        if ((my_pos == Child::cLU) || (my_pos == Child::cRU)) {
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
    if (n == CornerNeigh::LD) {
        if (my_pos != Child::cRU) {
            return ret;
        }
    }
    if (n == CornerNeigh::LU) {
        if (my_pos != Child::cRD) {
            return ret;
        }
    }
    if (n == CornerNeigh::LD) {
        if (my_pos != Child::cRU) {
            return ret;
        }
    }
    if (n == CornerNeigh::LD) {
        if (my_pos != Child::cRU) {
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
    if (n == Neigh::DOWN) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cLU));
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cRU));
        return halfSizeNeighs;
    }
    if (n == Neigh::UP) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cLD));
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cRD));
        return halfSizeNeighs;
    }
    if (n == Neigh::LEFT) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cRU));
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cRD));
        return halfSizeNeighs;
    }
    if (n == Neigh::RIGHT) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cLU));
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cLD));
        return halfSizeNeighs;
    }
    return halfSizeNeighs;
}

vector<TreeIndex>& TreeIndex::get_halfsize_possible_face_neighbours_optimized(vector<TreeIndex>& buf, Neigh n) {
    TreeIndex fullSizeNeigh = get_face_neighbor(n);
    if (n == Neigh::DOWN) {
        buf.push_back(fullSizeNeigh.get_child(Child::cLU));
        buf.push_back(fullSizeNeigh.get_child(Child::cRU));
        return buf;
    }
    if (n == Neigh::UP) {
        buf.push_back(fullSizeNeigh.get_child(Child::cLD));
        buf.push_back(fullSizeNeigh.get_child(Child::cRD));
        return buf;
    }
    if (n == Neigh::LEFT) {
        buf.push_back(fullSizeNeigh.get_child(Child::cRU));
        buf.push_back(fullSizeNeigh.get_child(Child::cRD));
        return buf;
    }
    if (n == Neigh::RIGHT) {
        buf.push_back(fullSizeNeigh.get_child(Child::cLU));
        buf.push_back(fullSizeNeigh.get_child(Child::cLD));
        return buf;
    }

    return buf;
}

vector<TreeIndex> TreeIndex::get_halfsize_possible_corner_neighbours(CornerNeigh n) {
    TreeIndex fullSizeNeigh = get_corner_neighbor(n);
    vector<TreeIndex> halfSizeNeighs;
    if (n == CornerNeigh::LU) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cRD));
        return halfSizeNeighs;
    }
    if (n == CornerNeigh::RU) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cLD));
        return halfSizeNeighs;
    }
    if (n == CornerNeigh::LD) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cRU));
        return halfSizeNeighs;
    }
    if (n == CornerNeigh::RD) {
        halfSizeNeighs.push_back(fullSizeNeigh.get_child(Child::cLU));
        return halfSizeNeighs;
    }
    return halfSizeNeighs;
}

vector<TreeIndex> TreeIndex::get_all_halfsize_possible_neighs() {
    vector<TreeIndex> res;
    res.reserve( 8 );

    if (!is_left_border()) {
        res = get_halfsize_possible_face_neighbours_optimized(res, Neigh::LEFT);
    }
    if (!is_right_border()) {
        res = get_halfsize_possible_face_neighbours_optimized(res, Neigh::RIGHT);
    }
    if (!is_upper_border()) {
        res = get_halfsize_possible_face_neighbours_optimized(res, Neigh::UP);
    }
    if (!is_down_border()) {
        res = get_halfsize_possible_face_neighbours_optimized(res, Neigh::DOWN);
    }

    return res;
}

vector<TreeIndex>& TreeIndex::get_all_halfsize_possible_neighs_optimized(vector<TreeIndex>& buf) {

    if (!is_left_border()) {
        buf = get_halfsize_possible_face_neighbours_optimized(buf, Neigh::LEFT);
    }
    if (!is_right_border()) {
        buf = get_halfsize_possible_face_neighbours_optimized(buf, Neigh::RIGHT);
    }
    if (!is_upper_border()) {
        buf = get_halfsize_possible_face_neighbours_optimized(buf, Neigh::UP);
    }
    if (!is_down_border()) {
        buf = get_halfsize_possible_face_neighbours_optimized(buf, Neigh::DOWN);
    }

    return buf;
}

vector<TreeIndex> TreeIndex::get_all_samesize_possible_neighs() {
    vector<TreeIndex> res;
    res.reserve( 4 );

    if (!is_left_border()) {
        res.push_back(get_face_neighbor(Neigh::LEFT));
    }
    if (!is_right_border()) {
        res.push_back(get_face_neighbor(Neigh::RIGHT));
    }
    if (!is_down_border()) {
        res.push_back(get_face_neighbor(Neigh::DOWN));
    }
    if (!is_upper_border()) {
        res.push_back(get_face_neighbor(Neigh::UP));
    }

    return res;
}

vector<TreeIndex>& TreeIndex::get_all_samesize_possible_neighs_optimized(vector<TreeIndex>& buf) {

    if (!is_left_border()) {
        buf.push_back(get_face_neighbor(Neigh::LEFT));
    }
    if (!is_right_border()) {
        buf.push_back(get_face_neighbor(Neigh::RIGHT));
    }
    if (!is_down_border()) {
        buf.push_back(get_face_neighbor(Neigh::DOWN));
    }
    if (!is_upper_border()) {
        buf.push_back(get_face_neighbor(Neigh::UP));
    }

    return buf;
}

vector<TreeIndex> TreeIndex::get_all_larger_possible_neighs() {
    vector<TreeIndex> res;
    res.reserve(4);

    if (!is_left_border()) {
        res = get_larger_possible_face_neighbour_optimized(res, Neigh::LEFT);
    }
    if (!is_right_border()) {
        res = get_larger_possible_face_neighbour_optimized(res, Neigh::RIGHT);
    }
    if (!is_upper_border()) {
        res = get_larger_possible_face_neighbour_optimized(res, Neigh::UP);
    }
    if (!is_down_border()) {
        res = get_larger_possible_face_neighbour_optimized(res, Neigh::DOWN);
    }

    return res;
}

vector<TreeIndex>& TreeIndex::get_all_larger_possible_neighs_optimized(vector<TreeIndex>& buf) {

    if (!is_left_border()) {
        buf = get_larger_possible_face_neighbour_optimized(buf, Neigh::LEFT);
    }
    if (!is_right_border()) {
        buf = get_larger_possible_face_neighbour_optimized(buf, Neigh::RIGHT);
    }
    if (!is_upper_border()) {
        buf = get_larger_possible_face_neighbour_optimized(buf, Neigh::UP);
    }
    if (!is_down_border()) {
        buf = get_larger_possible_face_neighbour_optimized(buf, Neigh::DOWN);
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
    for (TreeIndex c: all_neighs) {
        // cout << "(" << c.lvl << "," << c.i << "," << c.j << "), ";
        all_neighs_ids.push_back(c.get_global_number());
    }
    // cout << " }" << endl;

    sort(all_neighs_ids.begin(), all_neighs_ids.end());
    auto last = std::unique(all_neighs_ids.begin(), all_neighs_ids.end());
    all_neighs_ids.erase(last, all_neighs_ids.end());

    return all_neighs_ids;
}

bool TreeIndex::is_left_border() {
    TreeIndex c = get_face_neighbor(Neigh::LEFT);
    return (c.j < 0);
}

bool TreeIndex::is_right_border() {
    TreeIndex c = get_face_neighbor(Neigh::RIGHT);
    return ((c.j & (-1 << max_lvl)) != 0);
}

bool TreeIndex::is_upper_border() {
    TreeIndex c = get_face_neighbor(Neigh::UP);
    return ((c.i & (-1 << max_lvl)) != 0);
}

bool TreeIndex::is_down_border() {
    TreeIndex c = get_face_neighbor(Neigh::DOWN);
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


void TreeIndex::get_corner_coords(double *x, double *y) {
//    double lvl_dx = min_dx * pow(2, max_lvl - lvl);
    *x = min_dx * i;
    *y = min_dx * j;
}


/* * * * * * * * РАБОТА С ЯЧЕЙКАМИ И ДЕРЕВОМ * * * * * * * */

double dist(double x1, double y1, double x2, double y2) {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

void Cell::get_spacial_coords(double *x, double *y) {
    double lvl_dx = min_dx * pow(2, max_lvl - lvl);
    *x = min_dx * i + lvl_dx/2; 
    *y = min_dx * j + lvl_dx/2; 
}

void Cell::get_border_cond(char *cond_type, double (**cond_func)(double, double, double)) {
    if (is_left_border()) {
        Area::get_border_cond(Area::Border::LEFT, cond_type, cond_func);
        return;
    }

    if (is_right_border()) {
        Area::get_border_cond(Area::Border::RIGHT, cond_type, cond_func);
        return;
    }

    if (is_down_border()) {
        Area::get_border_cond(Area::Border::DOWN, cond_type, cond_func);
        return;
    }

    if (is_upper_border()) {
        Area::get_border_cond(Area::Border::UP, cond_type, cond_func);
        return;
    }

    *cond_type = -1;
}



int LinearTree::FindCell(GlobalNumber_t target, Cell *cell) {
    // cout << "FindCell(" << target << ")\n";

    // binary search
    int left = 0;
	int right = cells.size();
	while (left < right) {
		int midi = (right + left) / 2;
		Cell mid = cells[midi];
        GlobalNumber_t val = mid.get_global_number();
		if (val == target) {
			// return midi;
            if (cell != nullptr) {
                *cell = mid;
            }
            return midi;
		}
        if (val < target) {
			left = midi + 1;
		} else {
			right = midi;
		}
	}
	return -1;
}
void LinearTree::Write(string filename) {
    vector<char> buf = GenWriteStruct();
    int len = buf.size();

    std::ofstream fout(filename, std::ios::out | std::ios::binary);
    fout.write(&buf[0], len);
    fout.close();
}

void LinearTree::WriteOffsets(string filename, int n_of_procs) {
    int one_sz = 3 * sizeof(int) + sizeof(double);
    
    vector<int> offsets;
    int sum = 0;
    for (int i = 0; i < n_of_procs-1; i++) {
        offsets.push_back(sum * one_sz);
        int len = cells.size() / n_of_procs;
        offsets.push_back(len * one_sz);
        sum += len;
    }
    offsets.push_back(sum * one_sz);
    offsets.push_back((cells.size() - sum) * one_sz);

    cout << "OFFSETS={ ";
    for (int i = 0; i < offsets.size(); i++) {
        cout << offsets[i] << " ";
    }
    cout << "}" << endl;

    std::ofstream fout(filename, std::ios::out | std::ios::binary);
    fout.write((char *)&offsets[0], offsets.size() * sizeof(int));
    fout.close();
}

void LinearTree::GenFromWriteStruct(vector<char>& buf) {
    max_present_lvl = base_lvl;

    int one_sz = 3 * sizeof(int) + sizeof(double);

    char *p = &buf[0];
    int pos = 0;

    int lvl_offset = 0;
    int i_offset = sizeof(int);
    int j_offset = 2 * sizeof(int);
    int temp_offset = 3 * sizeof(int);

    while (pos < buf.size()) {
        Cell c;
        c.lvl  = * ((int *)(&p[pos+lvl_offset]));
        c.i    = * ((int *)(&p[pos+i_offset]));
        c.j    = * ((int *)(&p[pos+j_offset]));
        c.temp[0] = * ((double *)(&(p[pos+temp_offset])));
        c.refine_mark = 0;
        cells.push_back(c);

        // std::cout << "cell pushed pos=" << pos << endl; 
        if (c.lvl > max_present_lvl) {
            max_present_lvl = c.lvl;
        }

        pos += one_sz;
    }

    cout << "first cell = Cell(" << cells[0].lvl << ", " << cells[0].i << "," << cells[0].j << ", " << cells[0].temp[0] << ")";
    cout << "last cell = Cell(" << cells[cells.size()-1].lvl << ", " << cells[cells.size()-1].i << "," << cells[cells.size()-1].j << ", " << cells[cells.size()-1].temp[0] << ")";
}

double get_lvl_dx(int lvl) {
    return  min_dx * pow(2, max_lvl - lvl);
}



/* * * * * * * * РАБОТА БЛОЧНЫМ ДЕРЕВОМ * * * * * * * */


BlockOfCells::BlockOfCells(int _cells_lvl, int _blk_lvl, int _sz, GlobalNumber_t _i):
        cells_lvl(_cells_lvl), i(_i), idx(_blk_lvl, _i), sz(_sz) {}

void BlockOfCells::CreateCells(double (*Temp_func)(double, double)) {

    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz; j++) {

            double x, y;
            get_spacial_coords(i, j, &x, &y);

            SimpleCell c(Temp_func(x, y));
            cells.push_back(c);
        }
    }
}

void BlockOfCells::get_spacial_coords(int i, int j, double *x, double *y) {
    double xx, yy;
    idx.get_corner_coords(&xx, &yy);

    double lvl_dx = min_dx * pow(2, max_lvl - cells_lvl);
    xx += i * lvl_dx + lvl_dx/2;
    yy += j *lvl_dx + lvl_dx/2;

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
                cout << int(cells[i].refine_mark)  << " ";
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

    BlockOfCells bc00 = BlockOfCells(cells_lvl, idx.lvl+1, sz/2, idx.get_child(Child::cLD).get_global_number());
    bc00.CreateCells(Temp_func);
    if (refine_marks[0]) {
        bc00.refine_mark = 1;
    }
    children.push_back(bc00);

    BlockOfCells bc01 = BlockOfCells(cells_lvl, idx.lvl+1, sz/2, idx.get_child(Child::cRD).get_global_number());
    bc01.CreateCells(Temp_func);
    if (refine_marks[1]) {
        bc01.refine_mark = 1;
    }
    children.push_back(bc01);

    BlockOfCells bc10 = BlockOfCells(cells_lvl, idx.lvl+1, sz/2, idx.get_child(Child::cLU).get_global_number());
    bc10.CreateCells(Temp_func);
    if (refine_marks[2]) {
        bc10.refine_mark = 1;
    }
    children.push_back(bc10);

    BlockOfCells bc11 = BlockOfCells(cells_lvl, idx.lvl+1, sz/2, idx.get_child(Child::cRU).get_global_number());
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


BlockedLinearTree::BlockedLinearTree(double (*Temp_func)(double, double)) {
    max_present_lvl = base_lvl;
    max_present_blk_lvl = base_blk_lvl;

    int global_i_start = 0;
    int global_i_stop  = 1 << (max_blk_lvl*2);
    int global_i_step  = 1 << (2*max_blk_lvl - 2*base_blk_lvl);


    for (int i = global_i_start; i < global_i_stop; i += global_i_step) {
        BlockOfCells bc = BlockOfCells(base_lvl, base_blk_lvl, base_blk_sz, i);
        bc.CreateCells(Temp_func);
        blocks.push_back(bc);
    }
}

int BlockedLinearTree::MarkToRefine() {
    if (max_present_lvl == max_lvl) {
        return 0;
    }

    int n_of_block_marks = 0;
    for (int i = 0; i < blocks.size(); i++) {
        n_of_block_marks += blocks[i].MarkToRefine();
    }

    return (n_of_block_marks > 0);
}

void BlockedLinearTree::RefineBlocks() {

    int new_max_present_blk_lvl = max_present_blk_lvl;
    for (int i = 0; i < blocks.size();) {
        if (blocks[i].refine_mark) {
            if (blocks[i].idx.lvl < max_blk_lvl) {
                vector<BlockOfCells> new_blocks = blocks[i].Split(&Area::T0);
                blocks.erase(blocks.begin()+i);
                blocks.insert(blocks.begin()+(i), new_blocks.begin(), new_blocks.end());

                if (blocks[i].idx.lvl + 1 > new_max_present_blk_lvl) {
                    new_max_present_blk_lvl = blocks[i].idx.lvl + 1;
                }
                i+= 4;

            } else {
                i++;
            }
        }
    }

    if (new_max_present_blk_lvl > max_present_blk_lvl) {
        max_present_blk_lvl = new_max_present_blk_lvl;
    }
}

void BlockedLinearTree::RefineCells() {
    int new_max_present_lvl = max_present_lvl;

    for (int i = 0; i < blocks.size(); i++) {
        if (blocks[i].refine_mark) {
            blocks[i].RefineCells();
            if (blocks[i].cells_lvl > new_max_present_lvl) {
                new_max_present_lvl = blocks[i].cells_lvl;
            }
        }
    }

    max_present_lvl = new_max_present_lvl;
}

void BlockedLinearTree::Decompose(int n_procs) {

    int sum_weight = 0;
    for (int i = 0; i < blocks.size(); i++) {
        sum_weight += blocks[i].GetNOfCells();
    }

    proc_blocks = vector<int>(blocks.size(), 0);

    int optimal_weight = sum_weight / n_procs;
    int cur_weight = 0;
    int proc_i = 0;
    for (int i = 0; i < blocks.size(); i++) {
        proc_blocks[i] = proc_i;
        cur_weight += blocks[i].GetNOfCells();
        if (cur_weight >= optimal_weight) {
            proc_i++;
            cur_weight = 0;
        }
    }
}

void BlockedLinearTree::Write(string filename) {
    vector<char> buf = GenWriteStruct();
    int len = buf.size();

    std::ofstream fout(filename, std::ios::out | std::ios::binary);
    fout.write(&buf[0], len);
    fout.close();
}

vector<char> BlockedLinearTree::GenWriteStruct() {
    vector<char> buf;
    std::cout << " gen structs for write start\n";
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
            GlobalNumber_t glob_idx = get_glob_idx(blocks[blk_i].i, k, lvl);
            split_ints(glob_idx, &tree_i, &tree_j);
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

            int proc_n = proc_blocks[blk_i];
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
    std::cout << " gen structs for write finished\n";
    return buf;
}

GlobalNumber_t get_glob_idx(GlobalNumber_t blk_i, GlobalNumber_t cell_i, int cell_lvl) {
    GlobalNumber_t res = blk_i << 2*(max_lvl - max_blk_lvl);
    res = res | (cell_i << 2*(max_lvl - (max_blk_lvl + cell_lvl)));
    return res;
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

