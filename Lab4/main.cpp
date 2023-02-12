#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits.h>
#include <algorithm>
// #include "io.hpp"
using namespace std;

int N = 0, height = 21, snc = 10, nspan = 2;
vector<int> T, B;
vector<int> track;
vector<int> gm;
vector<int> last;
vector<vector<int>> Y;
vector<vector<int>> V;
vector<vector<vector<int>>> H;

struct pattern {
    vector<int> wires = {};
    int track = 0;
    int gap_to_side = 0;
    int jog = 0;
};

void parse(string& infile)
{
    ifstream fin(infile);
    string line, tmp, tmp_;
    int n;

    getline(fin, line);
    stringstream ssT(line);
    while (getline(ssT, tmp, ' ') && tmp != "\r")
    {
        stringstream ss_(tmp);
        while (getline(ss_, tmp_, '\t'))
        {
            n = stoi(tmp_);
            T.push_back(n);

            if (n > N) N = n;
        }
    }

    getline(fin, line);
    stringstream ssB(line);
    while (getline(ssB, tmp, ' ') && tmp != "\r")
    {
        stringstream ss_(tmp);
        while (getline(ss_, tmp_, '\t'))
        {
            n = stoi(tmp_);
            B.push_back(n);

            if (n > N) N = n;
        }
    }

    fin.close();
}

void initialize()
{
    last.resize(N+1, -1);
    Y.resize(N+1, {});
    track.resize(height, -1);
    gm.resize(height, -1);
    V.resize(height, vector<int>(T.size(), -1));
    H.resize(height);

    for (int i=T.size()-1; i>=0; i--)
    {
        if (last[T[i]] == -1) last[T[i]] = i;
        if (last[B[i]] == -1) last[B[i]] = i;
    }

    // int t = T[0], b = B[0];
    // bool tt = false, tb = false, bt = false, bb = false;
    // for (int i=1; i<snc; i++)
    // {
    //     if (T[i] == t) tt = true;
    //     if (T[i] == b) bt = true;
    //     if (B[i] == t) tb = true;
    //     if (B[i] == b) bb = true;
    // }

    // int ts = 1, bs = 1;
    // if (tt && !tb) ts = 2;
    // if (!tt && tb) ts = 0;
    // if (bt && !bb) bs = 2;
    // if (!bt && bb) bs = 0;

    // int idx_t = -1, idx_b = -1;
    // if (t != 0)
    // {
    //     if (ts == 2 || (b != 0 && bs == 2 && bs > ts)) {idx_t = height-2;}
    //     else if (ts == 1 || (ts == 0 && bs == 2)) {idx_t = height/2+1;}
    //     else if (ts == 0 && b != 0) {idx_t = 2;}
    //     else {idx_t = 1;}

    //     Y[t].push_back(idx_t);
    //     track[Y[t].back()] = t;
    //     H[Y[t].back()].emplace_back(vector<int>({0, t, 0}));
    // }

    // if (b != 0)
    // {
    //     if (bs == 0 || (t != 0 && ts == 0 && ts < bs)) {idx_b = 1;}
    //     else if (bs == 1 || (ts == 0 && bs == 2)) {idx_b = height/2;}
    //     else if (bs == 2 && t != 0) {idx_b = height-3;}
    //     else {idx_b = height-2;}

    //     Y[b].push_back(idx_b);
    //     track[Y[b].back()] = b;
    //     H[Y[b].back()].emplace_back(vector<int>({0, b, 0}));
    // }
}

void make_feasible_top_and_bottom_connections(int col)
{
    int t = T[col], b = B[col];

    // top and bottom signal are the same (Ti == Bi)
    if (t == b && t != 0)
    {
        for (int i=0; i<height; i++) 
        {
            if (track[i] == t) track[i] = -1;
            gm[i] = t;
        }
        Y[t].clear();
        return;
    }

    // search for target track for Ti and Bi
    int t_ = INT_MAX, b_ = -INT_MAX;
    if (t != 0)
    {
        for (int i=height-2; i>0; i--)
        {
            if (track[i] == t || track[i] == -1)
            {
                t_ = i;
                break;
            }
        }
    }

    if (b != 0)
    {
        for (int i=1; i<min(t_, height-1); i++)
        {
            if (track[i] == b || track[i] == -1)
            {
                b_ = i;
                break;
            }
        }
    }

    if (b_ > t_)
    {
        // top and bottom vertical segments overlap
        if (height - t_ > b_)
        {
            // choose bottom vertical segment
            if (track[b_] == -1) 
            {
                Y[b].push_back(b_);
                track[b_] = b;
                H[b_].emplace_back(vector<int>({col, b, col}));
            }
            track[0] = -1;
            for (int i=0; i<=b_; i++) gm[i] = b;
        }
        else
        {
            // choose top vertical segment
            if (track[t_] == -1) 
            {
                Y[t].push_back(t_);
                track[t_] = t;
                H[t_].emplace_back(vector<int>({col, t, col}));
            }
            track[height-1] = -1;
            for (int i=t_; i<=height-1; i++) gm[i] = t;
        }
    }
    else
    {
        if (t_ != INT_MAX)
        {
            if (track[t_] == -1) 
            {
                Y[t].push_back(t_);
                track[t_] = t;
                H[t_].emplace_back(vector<int>({col, t, col}));
            }
            else if (Y[t].size() == 1 && last[t] <= col)
            {
                track[t_] = -1;
                Y[t].clear();
            }
            track[height-1] = -1;
            for (int i=t_; i<=height-1; i++) gm[i] = t;
        }
        
        if (b_ != -INT_MAX)
        {
            if (track[b_] == -1) 
            {
                Y[b].push_back(b_);
                track[b_] = b;
                H[b_].emplace_back(vector<int>({col, b, col}));
            }
            else if (Y[b].size() == 1 && last[b] <= col)
            {
                track[b_] = -1;
                Y[b].clear();
            }
            track[0] = -1;
            for (int i=0; i<=b_; i++) gm[i] = b;
        }
    }
}

bool cmp(pattern& tmp, pattern& best)
{
    if (tmp.track < best.track) return false;
    if (tmp.track == best.track)
    {
        if (tmp.gap_to_side < best.gap_to_side || 
            (tmp.gap_to_side == best.gap_to_side && tmp.jog > best.jog))
            return true;
        else
            return false;
    }
    return true;
}

void search_patterns(vector<vector<int>>& seg, pattern& tmp, pattern& best, int idx)
{
    if (idx == seg.size()) 
    {
        if (tmp.wires.size() > 0)
            tmp.gap_to_side = (seg[tmp.wires[0]][0]) + (height - seg[tmp.wires.back()][1]);
        if (cmp(tmp, best)) best = tmp;
        return;
    }

    // the pattern does not pick this vertical segment
    search_patterns(seg, tmp, best, idx+1);

    // the pattern picks this vertical segment
    if (tmp.wires.size() == 0 ||
        seg[tmp.wires.back()][1] < seg[idx][0] ||
        (seg[tmp.wires.back()][1] == seg[idx][0] && track[seg[tmp.wires.back()][1]] == track[seg[idx][0]]))
    {
        // check if the segment would free two tracks
        int tr = 1;
        if (seg[idx][3] == 1 && 
            !(tmp.wires.size() != 0 && seg[tmp.wires.back()][1] == seg[idx][0]))
                tr = 2;
        
        tmp.track += tr;
        tmp.wires.push_back(idx);
        tmp.jog += seg[idx][2];
        search_patterns(seg, tmp, best, idx+1);
        tmp.jog -= seg[idx][2];
        tmp.wires.pop_back();
        tmp.track -= tr;
    }
}

void free_track(int r)
{
    auto it_y = find(Y[track[r]].begin(), Y[track[r]].end(), r);
    if (it_y != Y[track[r]].end()) Y[track[r]].erase(it_y);
    track[r] = -1;
}

void free_tracks_as_possible(int col)
{
    vector<vector<int>> seg;
    pattern tmp, best;
    int bot = B[col] == 0 ? 1 : 0, top = T[col] == 0 ? height-1 : height;

    // compute all feasible vertical segments
    for (int i=bot; i<top; i++)
    {
        // track i is already occupied by a segment
        if (gm[i] != -1)
        {
            int idx = gm[i++];
            while (i < height && gm[i] == idx) {i++;}
            i--;
        }

        if (track[i] == -1 || Y[track[i]].size() < 2) continue;

        for (int j=i+1; j<height; j++)
        {
            if (track[j] == track[i])
            {
                int end = last[track[i]] > col+nspan ? 0 : 1;
                seg.emplace_back(vector<int>({i, j, j-i, end}));
                break;
            }
            else if (gm[j] != -1)
                break;
        }
    }

    // search all feasible patterns
    search_patterns(seg, tmp, best, 0);

    // update tracks
    int down, up;
    vector<int>::iterator it_y;
    for (auto it=best.wires.begin(); it!=best.wires.end(); it++)
    {
        down = seg[*it][0];
        up = seg[*it][1];
        for (int i=down; i<=up; i++) gm[i] = track[up];
      
        // the segment free two tracks
        if (last[track[up]] <= col && Y[track[up]].size() == 2) 
        {
            free_track(up);
            free_track(down);
        }
        else
        {
            int n = track[up], flag = 0;
            for (int i=col+1; i<T.size(); i++)
            {
                if (T[i] == n) 
                {
                    free_track(down);
                    flag = 1;
                    break;
                }
                else if (B[i] == n)
                {
                    free_track(up);
                    flag = 1;
                    break;
                }
            }

            if (flag == 0) free_track(down);
        }
    }
}

bool near_signal_exist(int tmp, int col, vector<int>& sig)
{
    int to = min((int)T.size(), col+nspan);
    for (int i=col+1; i<to; i++)
        if (sig[i] == tmp) 
            return true;
    return false;
}

void add_jogs_to_reduce_range_of_split_nets(int col)
{
    for (int i=0; i<Y.size(); i++)
    {
        if (Y[i].size() < 2) continue;

        int min = Y[i][0], max = Y[i][0];
        for (int j=1; j<Y[i].size(); j++)
        {
            if (Y[i][j] < min) min = Y[i][j];
            else if (Y[i][j] > max) max = Y[i][j];
        }

        int to = min+1;
        if (!near_signal_exist(i, col, B))
        {
            while (to < height && (gm[to] == -1 || gm[to] == i)) {to++;}
            while (--to > min && (track[to] != -1 || gm[to] == i)) {}
            if (to != min)
            {
                for (int j=min; j<=to; j++) {gm[j] = i;}
                H[to].emplace_back(vector<int>({col, i, col}));
                free_track(min);
                Y[i].push_back(to);
                track[to] = i;
            }
        }
        
        to = max-1;
        while (to >=0 && (gm[to] == -1 || gm[to] == i)) {to--;}
        while (++to < max && (track[to] != -1 || gm[to] == i)) {}
        if (to != max)
        {
            for (int j=max; j>=to; j--) {gm[j] = i;}
            H[to].emplace_back(vector<int>({col, i, col}));
            free_track(max);
            Y[i].push_back(to);
            track[to] = i;
        }
    }
}

void add_jogs_to_raise_rising_and_falling_nets(int col)
{
    vector<int> dist(height, 0);
    vector<int> order;
    int lookTo = min(col+snc, (int)T.size());
    int bot = B[col] == 0 ? 1 : 0, top = T[col] == 0 ? height-1 : height;

    // calcutate distance
    for (int i=bot; i<top; i++)
    {
        // not unsplit net
        if (track[i] == -1 || Y[track[i]].size() > 1) continue;
        
        // check if the net is rising or falling
        bool t = false, b = false;
        int n = track[i];
        for (int j=col+1; j<lookTo && !t && !b; j++)
        {
            if (T[j] == n) t = true;
            if (B[j] == n) b = true;
        }

        // set distance
        int to = i;
        if (t && !b) 
        {
            while (++to < height && gm[to] == -1) {}
            while (--to > i && track[to] != -1) {}
            dist[i] = to - i;
        }
        else if (!t && b) 
        {
            while (--to >= 0 && gm[to] == -1) {}
            while (++to < i && track[to] != -1) {}
            dist[i] = to - i;
        }
    }

    // order of decreasing distance
    order.reserve(height);
    for (int i=0; i<height; i++) order.push_back(i);
    std::sort(order.begin(), order.end(), [&](int a, int b){
        return abs(dist[a]) > abs(dist[b]);
    });

    // find jog
    for (auto it=order.begin(); it!=order.end(); it++)
    {
        if (gm[*it] != -1 && gm[*it] != track[*it]) continue;

        int to = *it;
        if (dist[*it] > 0)
        {
            // upward jog
            while (++to < height && gm[to] == -1) {}
            while (--to > *it && track[to] != -1) {}
        }
        else if (dist[*it] < 0)
        {
            // downward jog
            while (--to >= 0 && gm[to] == -1) {}
            while (++to < *it && track[to] != -1) {}
        }
        else
        {
            break;
        }

        if (*it != to)
        {
            Y[track[*it]].push_back(to);
            H[to].emplace_back(vector<int>({col, track[*it], col}));
            track[to] = track[*it];
            free_track(*it);

            int f = min(*it, to), t = max(*it, to);
            for (int i=f; i<=t; i++) gm[i] = track[to];
        }
    }    
}

void update_Y(int idx)
{
    for (int i=0; i<Y.size(); i++)
    {
        for (int j=0; j<Y[i].size(); j++)
            if (Y[i][j] >= idx)
                Y[i][j]++;
    }
}

void update_channel(int col)
{
    int idx = -2, from, cnt = 0;
    for (int i=0; i<height; i++)
    {
        if (gm[i] != -1)
        {
            idx = gm[i];
            from = i++;

            while (i < height && gm[i] == idx) {i++;}
            i--;

            V[from][col] = V[i][col] = idx;
        }
    }

    if (gm[height-1] == -1 && T[col] != 0)
    {
        for (int i=height-2; i>0; i--)
        {
            if (gm[i] != -1 || i==height/2+1)
            {
                V.insert(V.begin()+(i+1), vector<int>(T.size(), -1));
                H.insert(H.begin()+(i+1), vector<vector<int>>({}));
                V[height][col] = V[i+1][col] = T[col];
                H[i+1].emplace_back(vector<int>({col, T[col], col}));
                update_Y(i+1);
                Y[T[col]].push_back(i+1);
                track.insert(track.begin()+(i+1), T[col]);
                break;
            }
        }
        cnt++;
    }

    if (gm[0] == -1 && B[col] != 0)
    {
        for (int i=1; i<height-1; i++)
        {
            if (gm[i] != -1 || i==height/2)
            {
                V.insert(V.begin()+i, vector<int>(T.size(), -1));
                H.insert(H.begin()+i, vector<vector<int>>({}));
                V[0][col] = V[i][col] = B[col];
                H[i].emplace_back(vector<int>({col, B[col], col}));
                update_Y(i);
                Y[B[col]].push_back(i);
                track.insert(track.begin()+i, B[col]);
                break;
            }
        }
        cnt++;
    }

    if (T[col] == B[col] && T[col] != 0 && last[T[col]] > col)
    {
        int from_, to_ = -1;
        if (T[col] == T[last[T[col]]])
        {
            from_ = height+cnt-2;
            for (int i=from_; i>0; i--)
            {
                if (track[i] == -1)
                {
                    to_ = i;
                    break;
                }
            }
        }
        else
        {
            from_ = 1;
            for (int i=from_; i<height+cnt-1; i++)
            {
                if (track[i] == -1)
                {
                    to_ = i;
                    break;
                }
            }
        }

        if (to_ == -1)
        {
            to_ = T[col] == T[last[T[col]]] ? height-1+cnt : 1;
            update_Y(to_);
            V.insert(V.begin()+to_, vector<int>(T.size(), -1));
            H.insert(H.begin()+to_, vector<vector<int>>({}));
            track.insert(track.begin()+to_, T[col]);
            cnt++;
        }
        H[to_].emplace_back(vector<int>({col, T[col], col}));
        V[from][col] = V[to_][col] = T[col];
        Y[T[col]].push_back(to_);
        track[to_] = T[col];
    }

    height += cnt;
    gm.resize(height, -1);
    for (int i=0; i<height; i++) gm[i] = -1;
}

void extend_to_next_column(int col)
{
    for (int i=1; i<height-1; i++)
    {
        if (track[i] != -1)
        {
            if (Y[track[i]].size() <= 1 && last[track[i]] <= col)
                Y[track[i]].clear();
            else
                H[i].back()[2]++;
        }
    }
}

void output(string& outfile)
{
    vector<vector<string>> Net(N+1, vector<string>({}));

    // horizontal
    for (int i=1; i<height-1; i++)
    {
        for (auto it=H[i].begin(); it!=H[i].end(); it++)
        {
            if ((*it)[0] == (*it)[2]) continue;
            Net[(*it)[1]].emplace_back(".H " + to_string((*it)[0]) + " " + to_string(i) + " " + to_string((*it)[2]));
        }
    }

    // vertical
    for (int j=0; j<T.size(); j++)
    {
        int i = 0, idx = -2, from = 0;

        while (i < height)
        {
            if (V[i][j] == idx) 
            {
                Net[idx].emplace_back(".V " + to_string(j) + " " + to_string(from) + " " + to_string(i));
                from = i;
            }

            if (V[i][j] != -1 && V[i][j] != idx) 
            {
                idx = V[i][j];
                from = i;
            }

            i++;
        }
    }

    ofstream fout(outfile);

    for (int i=0; i<Net.size(); i++)
    {
        if (Net[i].size() == 0) continue;

        fout << ".begin " << i << endl;
        
        for (int j=0; j<Net[i].size(); j++)
            fout << Net[i][j] << endl;

        fout << ".end" << endl;
    }

    fout.close();
}

bool greedy_finish(int col)
{
    if (col < T.size()) return false;

    for (int i=0; i<Y.size(); i++)
    {
        if (!Y[i].empty())
        {
            T.push_back(0);
            B.push_back(0);
            for (int j=0; j<height; j++) V[j].push_back(-1);
            return false;
        }
    }

    return true;
}

void greedy()
{
    initialize();

    int i = 0;
    while (!greedy_finish(i))
    {

        if (i == 13)
            srand(13);
        track[height-1] = T[i];
        track[0] = B[i];

        make_feasible_top_and_bottom_connections(i);

        free_tracks_as_possible(i);

        add_jogs_to_reduce_range_of_split_nets(i);

        add_jogs_to_raise_rising_and_falling_nets(i);

        update_channel(i);

        extend_to_next_column(i);

        i++;


        string out = "output.txt", draw = "outdraw.txt";
        output(out);
        io::drawNets("Case/case2.txt", out.c_str(), draw.c_str());
        // io::drawNets("Case/Deutsch_difficult.txt", out.c_str(), draw.c_str());
        srand(0);
    }


    // string out = "output.txt", draw = "outdraw.txt";
    // output(out);
    // io::drawNets("Case/Deutsch_difficult.txt", out.c_str(), draw.c_str());
    // srand(0);
}


int main(int argc, char* argv[])
{
    string infile = argv[1], outfile = argv[2];

    parse(infile);

    greedy();

    output(outfile);
}