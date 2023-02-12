#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <random>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <list>
#include <utility>

#define ll long long

using namespace std;

int outline_width, outline_height, num_blocks, num_terminals, num_nets;
int root, root_s;
int tmp_w, tmp_h, tmp_w_s, tmp_h_s, tmp_w_f, tmp_h_f;
ll HPWL, HPWL_s, best_cost, AREA, AREA_s, HPWL_f, AREA_f, cost_f = LONG_LONG_MAX, runtime = 289;
double alpha, T, epsilon = 1, r = 0.85, k = 15, up_limit = 40;
clock_t T_s;
map<string, int> idx_map;
vector<string> idx2name;
vector<int> bx, by, bw, bh;
vector<int> p, lc, rc;
vector<int> bx_s, by_s, bw_s, bh_s;
vector<int> p_s, lc_s, rc_s;
vector<int> bx_f, by_f, bw_f, bh_f;
list<pair<int, int>> contour;
vector<vector<int>> Net;
vector<double> pert_prob = {0.33, 0.66};

random_device rd;
mt19937 g(rd());
uniform_real_distribution<double> gen01(0.0, 1.0);
uniform_int_distribution<> genblk;

void parse(string& in_block, string& in_net)
{
    // .block
    ifstream fin_b(in_block);
    string tmp_x, tmp_y, name;

    fin_b >> tmp_x >> outline_width >> outline_height;
    fin_b >> tmp_x >> num_blocks;
    fin_b >> tmp_x >> num_terminals;

    bx.resize(num_blocks + num_terminals + 1, -1);
    by.resize(num_blocks + num_terminals + 1, -1);
    bw.resize(num_blocks + 1, -1);
    bh.resize(num_blocks + 1, -1);
    p.resize(num_blocks + 1, -1);
    lc.resize(num_blocks + 1, -1);
    rc.resize(num_blocks + 1, -1);
    bx_s.resize(num_blocks + num_terminals + 1, -1);
    by_s.resize(num_blocks + num_terminals + 1, -1);
    bw_s.resize(num_blocks + 1, -1);
    bh_s.resize(num_blocks + 1, -1);
    p_s.resize(num_blocks + 1, -1);
    lc_s.resize(num_blocks + 1, -1);
    rc_s.resize(num_blocks + 1, -1);
    bx_f.resize(num_blocks + num_terminals + 1, -1);
    by_f.resize(num_blocks + num_terminals + 1, -1);
    bw_f.resize(num_blocks + 1, -1);
    bh_f.resize(num_blocks + 1, -1);
    idx2name.resize(num_blocks + 1);
    genblk = uniform_int_distribution<>(1, num_blocks);

    // parse blocks
    for (int i=1; i<=num_blocks; i++)
    {
        fin_b >> name;
        idx2name[i] = name;
        idx_map[name] = i;
        fin_b >> bw[i] >> bh[i];
    }

    // parse terminals
    for (int i=num_blocks+1; i<=num_blocks + num_terminals; i++)
    {
        fin_b >> name >> tmp_x;
        idx_map[name] = i;
        fin_b >> bx[i] >> by[i];
    }
                                                                       
    fin_b.close();


    // .nets
    ifstream fin_n(in_net);
    string tmp_blk;
    int deg;

    fin_n >> tmp_blk >> num_nets;
    Net.resize(num_nets);

    for (int i=0; i<num_nets; i++)
    {
        fin_n >> tmp_blk >> deg;
        Net[i].resize(deg);

        for (int j=0; j<deg; j++)
        {
            fin_n >> tmp_blk;
            Net[i][j] = idx_map[tmp_blk];
        }
    }

    fin_n.close();
}

ll calCost()
{
    AREA = tmp_w * tmp_h;
    HPWL = 0;

    for (auto it=Net.begin(); it!=Net.end(); it++)
    {
        int b_ = INT_MAX, t_ = 0, l_ = INT_MAX, r_ = 0;

        for (auto itt=it->begin(); itt!=it->end(); itt++)
        {
            int x = bx[*itt], y = by[*itt];
            if (*itt <= num_blocks)
            {
                x +=  bw[*itt]/2;
                y +=  bh[*itt]/2;
            }

            if (x < l_) l_ = x;
            if (x > r_) r_ = x;
            if (y < b_) b_ = y;
            if (y > t_) t_ = y;
        }

        HPWL += ((r_ - l_) + (t_ - b_));
    }


    ll cost = (ll)((alpha * AREA) + ((1 - alpha) * HPWL));
    if (alpha == 0 && num_nets == 0) cost = AREA;

    if (tmp_w > outline_width || tmp_h > outline_height)
    {
        double alpha_w = max(1, 1 + tmp_w - outline_width), alpha_h = max(1, 1 + tmp_h - outline_height);
        return (ll)(alpha_w * alpha_h * cost);
    }
    else
        return cost;
}

void rotate()
{
    int n = genblk(g);
    swap(bw[n], bh[n]);
}

void pc_swap(int n1, int n2)
{
    int up = n1, down = n2;
    if (p[n1] == n2) swap(up, down);

    int l_ = lc[up], r_ = rc[up], p_ = p[up];

    if (l_  == down)
    {
        if (r_ != -1) p[r_] = l_;
        lc[p_] == up ? lc[p_] = l_ : rc[p_] = l_;
        if (lc[l_] != -1) p[lc[l_]] = up;
        if (rc[l_] != -1) p[rc[l_]] = up;
        p[l_] = p_;
        p[up] = l_;
        lc[up] = lc[l_];
        rc[up] = rc[l_];
        lc[l_] = up;
        rc[l_] = r_;
    }
    else
    {
        if (l_ != -1) p[l_] = r_;
        lc[p_] == up ? lc[p_] = r_ : rc[p_] = r_;
        if (lc[r_] != -1) p[lc[r_]] = up;
        if (rc[r_] != -1) p[rc[r_]] = up;
        p[r_] = p_;
        p[up] = r_;
        lc[up] = lc[r_];
        rc[up] = rc[r_];
        lc[r_] = l_;
        rc[r_] = up;
    }
    
    if (n1 == root) root = n2;
    else if (n2 == root) root = n1;
}

void reg_swap(int n1, int n2)
{
    int n1l = lc[n1], n1r = rc[n1], n1p = p[n1];
    int n2l = lc[n2], n2r = rc[n2], n2p = p[n2]; 
    bool left1 = lc[n1p] == n1 ? true : false;
    bool left2 = lc[n2p] == n2 ? true : false;

    if (n1l != -1) p[n1l] = n2;
    if (n1r != -1) p[n1r] = n2;
    if (n2l != -1) p[n2l] = n1;
    if (n2r != -1) p[n2r] = n1;

    left1 ? lc[n1p] = n2 : rc[n1p] = n2;
    left2 ? lc[n2p] = n1 : rc[n2p] = n1;

    swap(lc[n1], lc[n2]);
    swap(rc[n1], rc[n2]);
    swap(p[n1], p[n2]);

    if (root == n1) root = n2;
    else if (root == n2) root = n1;
}

void swap()
{
    if (num_blocks == 1) return;

    int n1 = genblk(g), n2 = genblk(g);
    while (n1 == n2) {n2 = genblk(g);}

    if (p[n1] == n2 || p[n2] == n1)
        pc_swap(n1, n2);
    else
        reg_swap(n1, n2);
}


void delete_b(int n)
{
    if (lc[n] != -1 && rc[n] == -1)
    {   
        if (n == root)
        {
            p[lc[n]] = 0;
            root = lc[n];
        }
        else
        {
            p[lc[n]] = p[n];
            lc[p[n]] == n ? lc[p[n]] = lc[n] : rc[p[n]] = lc[n];
        }
    }
    else if (lc[n] == -1 && rc[n] != -1)
    {
        if (n == root)
        {
            p[rc[n]] = 0;
            root = rc[n];
        }
        else
        {
            p[rc[n]] = p[n];
            rc[p[n]] == n ? rc[p[n]] = rc[n] : lc[p[n]] = rc[n];
        }
    }
    else if (lc[n] == -1 && rc[n] == -1)
    {
        lc[p[n]] == n ? lc[p[n]] = -1 : rc[p[n]] = -1;
    }
    else
    {
        while (true)
        {
            int l_ = lc[n], r_ = rc[n];

            if (l_ == -1 && r_ == -1) break;

            if ((gen01(g) < 0.5 && l_ != -1) || r_ == -1)
                pc_swap(n, l_);
            else
                pc_swap(n, r_);
        }

        lc[p[n]] == n ? lc[p[n]] = -1 : rc[p[n]] = -1;
    }

    p[n] = lc[n] = rc[n] = -1;
}

void random_insert(int blk, int to)
{
    int child;
    p[blk] = to;

    if (gen01(g) < 0.5)
    {
        child = lc[to];
        lc[to] = blk;
    }
    else
    {
        child = rc[to];
        rc[to] = blk;
    }

    gen01(g) < 0.5 ? lc[blk] = child : rc[blk] = child;
    if (child != -1) p[child] = blk;
}

void move()
{
    if (num_blocks == 1) return;

    int n1 = genblk(g), n2 = genblk(g);
    while (n1 == n2) {n2 = genblk(g);}

    delete_b(n1);
    random_insert(n1, n2);
}

int perturb()
{
    double rand = gen01(g);

    if (rand < pert_prob[0]) 
    {
        rotate(); 
        return 0;
    }
    else if (rand < pert_prob[1])
    {
        move(); 
        return 1;
    }
    else 
    {
        swap();
        return 2;
    }
}

void update_best()
{
    for (int i=0; i<=num_blocks; i++)
    {
        bx_s[i] = bx[i];
        by_s[i] = by[i];
        bw_s[i] = bw[i];
        bh_s[i] = bh[i];
        p_s[i] = p[i];
        lc_s[i] = lc[i];
        rc_s[i] = rc[i];
    }

    tmp_w_s = tmp_w;
    tmp_h_s = tmp_h;
    root_s = root;
    HPWL_s = HPWL;
    AREA_s = AREA;
}

void retrieve_best()
{
    for (int i=0; i<=num_blocks; i++)
    {
        bx[i] = bx_s[i];
        by[i] = by_s[i];
        bw[i] = bw_s[i];
        bh[i] = bh_s[i];
        p[i] = p_s[i];
        lc[i] = lc_s[i];
        rc[i] = rc_s[i];
    }

    tmp_w = tmp_w_s;
    tmp_h = tmp_h_s;
    root = root_s;
    HPWL = HPWL_s;
    AREA = AREA_s;
}

void update_to_final()
{
    for (int i=0; i<=num_blocks; i++)
    {
        bx_f[i] = bx_s[i];
        by_f[i] = by_s[i];
        bw_f[i] = bw_s[i];
        bh_f[i] = bh_s[i];
    }

    tmp_w_f = tmp_w_s;
    tmp_h_f = tmp_h_s;
    HPWL_f = HPWL_s;
    AREA_f = AREA_s;
    cost_f = best_cost;
}

void dfs(int tmp)
{
    bx[tmp] = lc[p[tmp]] == tmp ? bx[p[tmp]] + bw[p[tmp]] : bx[p[tmp]];
    int max_y = 0, r_end = bx[tmp] + bw[tmp];
    auto idx = contour.begin(), l_idx = idx;

    for (; idx->first < bx[tmp]; idx++) {}

    if (idx->first == bx[tmp] && idx->first != 0) 
        idx++;
    else
        contour.insert(idx, make_pair(bx[tmp], idx->second));

    l_idx = contour.insert(idx, make_pair(bx[tmp], -1));

    while(idx->first < r_end)
    {
        max_y = max(max_y, idx->second);
        contour.erase(idx++);
    }

    if (idx->first == r_end)
        contour.erase(idx++);
    else
        idx = contour.insert(idx, make_pair(r_end, idx->second));

    l_idx->second = max_y + bh[tmp];
    contour.insert(idx, make_pair(r_end, l_idx->second));

    by[tmp] = max_y;
    tmp_h = max(tmp_h, l_idx->second);
    if (r_end > tmp_w) tmp_w = r_end;

    if (lc[tmp] != -1) dfs(lc[tmp]);
    if (rc[tmp] != -1) dfs(rc[tmp]);
}

void packing()
{
    contour.clear();
    contour.push_back(make_pair(0, 0));
    contour.push_back(make_pair(INT_MAX, 0));
    tmp_w = tmp_h = 0;

    if (root > 0) dfs(root);
}

void update_pert_prob(int idx)
{
    switch(idx)
    {
        case 0:
            if (pert_prob[1] - pert_prob[0] >= 0.2 && pert_prob[1] <= 0.8)
            {
                pert_prob[0] += 0.01;
                pert_prob[1] += 0.005;
            }
        case 1:
            if (pert_prob[0] >= 0.2 && pert_prob[1] <= 0.8)
            {
                pert_prob[0] -= 0.01;
                pert_prob[1] += 0.05;
            }
        case 2:
            if (pert_prob[0] >= 0.2 && pert_prob[1] - pert_prob[0] >= 0.2)
            {
                pert_prob[0] -= 0.005;
                pert_prob[1] -= 0.01;
            }
    }
}

void SA()
{
    int MT = 1, uphill = 0, reject = 0;
    ll delta_cost = 0, tmp_cost;

    while (reject/MT < 0.95 && T > epsilon)
    {
        MT = uphill = reject = 0;

        while (uphill < up_limit && MT < k*num_blocks)
        {
            int idx = perturb();

            packing();
            
            tmp_cost = calCost();
            MT++;
            delta_cost = tmp_cost - best_cost;

            if (delta_cost <= 0 || gen01(g) < exp(-1 * delta_cost/T))
            {
                if (delta_cost > 0) 
                {
                    uphill++;
                }
                else 
                {
                    update_pert_prob(idx);
                }
                update_best();
                best_cost = tmp_cost;
            }
            else
            {
                reject++;
                retrieve_best();
            }
        }

        T *= r;
    }
}

void initBsTree()
{
    vector<int> blk_order;
    blk_order.resize(num_blocks);
    for (int i=1; i<=num_blocks; i++)  blk_order[i-1] = i;
    shuffle(blk_order.begin(), blk_order.end(), g);

    root = blk_order[0];
    p[root] = 0;
    bx[0] = by[0] = bw[0] = bh[0] = 0;
    int l_idx, r_idx, tmp_idx;

    for (int i=1; i<=num_blocks/2; i++)
    {
        tmp_idx = blk_order[i-1];

        if (2*i-1 < num_blocks)
        {
            l_idx =  blk_order[2*i-1];

            lc[tmp_idx] = l_idx;
            p[l_idx] = tmp_idx; 
        }
        else
        {
            lc[tmp_idx] = -1;
        }

        if (2*i < num_blocks)
        {
            r_idx = blk_order[2*i];
            
            rc[tmp_idx] = r_idx;
            p[r_idx] = tmp_idx;
        }
        else
        {
            rc[tmp_idx] = -1;
        }
    }

    for (int i=num_blocks/2+1; i<=num_blocks; i++)
    {
        tmp_idx = blk_order[i-1];
        lc[tmp_idx] = -1;
        rc[tmp_idx] = -1;
    }

    packing();
    update_best();
    best_cost = calCost();
    HPWL_s = HPWL;
    AREA_s = AREA;
    T = best_cost / 2;
    pert_prob = {0.3, 0.6};
}

void enable_outline_limit()
{
    ll tmp_cost, delta_cost;
    T = best_cost / 4;

    while(tmp_w > outline_width || tmp_h > outline_height)
    {
        int MT = 0, uphill = 0;

        while (uphill < up_limit && MT < k*num_blocks)
        {
            perturb();

            packing();
            
            tmp_cost = calCost();
            MT++;
            delta_cost = tmp_cost - best_cost;

            if (delta_cost <= 0 || gen01(g) < exp(-1 * delta_cost/T))
            {
                if (delta_cost > 0) uphill++;
                update_best();
                best_cost = tmp_cost;
            }
            else
                retrieve_best();

            if (tmp_w <= outline_width && tmp_h <= outline_height) break;
        }

        T *= r;
        if (T < 1) T = best_cost / 2;
    }

    if (alpha == 0 && num_nets == 0) best_cost = 0;
}

void output_graph_print()
{
    ofstream fout("draw.txt");

    fout << num_blocks << endl;
    fout << tmp_w << " " << tmp_h << endl;

    for (int i=1; i<=num_blocks; i++)
    {
        fout << i << " "
            << bx[i] << " "
            << by[i] << " "
            << bw[i] << " "
            << bh[i] << endl;
    }

    fout.close();
}

void output(string output)
{
    ofstream fout(output);

    fout << cost_f << endl;
    fout << HPWL_f << endl;
    fout << AREA_f << endl;
    fout << tmp_w_f << " " << tmp_h_f << endl;
    fout << ((double)(clock() - T_s)) / CLOCKS_PER_SEC << endl;

    for (int i=1; i<=num_blocks; i++)
    {
        fout << idx2name[i] << " "
            << bx_f[i] << " " << by_f[i] << " "
            << bx_f[i] + bw_f[i] << " " << by_f[i] + bh_f[i] << endl;
    }

    fout.close();
}

int main(int argc, char *argv[])
{
    T_s = clock();

    alpha = stod(argv[1]);
    string in_block = argv[2], in_net = argv[3], outfile = argv[4];
    int max_time = 0;
    clock_t pre_time = clock();

    parse(in_block, in_net);

    do
    {
        initBsTree();

        SA();

        enable_outline_limit();

        if (best_cost  < cost_f) update_to_final();

        max_time = max(max_time, (int)((clock() - pre_time) / CLOCKS_PER_SEC));
        pre_time = clock();
    } 
    while((int)((clock() - T_s) / CLOCKS_PER_SEC) + max_time < runtime);

    output(outfile);

    return 0;
}
