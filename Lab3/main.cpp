#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <time.h>
#include <future>
#include "abacus/abacus.hpp"
#include "FM/FM.hpp"
using namespace std;

int W, H, rh, rn;
int terNum, cellNum, virTerNum;
vector<string> idx2name;
vector<int> lx, ly, bw, bh, lay;
vector<int> lx0, ly0, bw0, bh0;
vector<int> lx1, ly1, bw1, bh1;
vector<int> ox, oy, olay;
vector<int> d0map, d1map;
vector<vector<pair<int, int>>> ter;


void parse(string& infile)
{
    ifstream fin(infile);
    if (!fin)
    {
        cout << "File could not be opened";
        fin.close();
        return;
    }

    string trash;

    fin >> trash >> W >> H;
    fin >> trash >> rh >> rn;
    
    fin >> trash >> terNum;
    lx.resize(terNum);
    ly.resize(terNum);
    bw.resize(terNum);
    bh.resize(terNum);
    idx2name.resize(terNum);
    for (int i=0; i<terNum; i++)
        fin >> idx2name[i] >> lx[i] >> ly[i] >> bw[i] >> bh[i];

    fin >> trash >> cellNum;
    lx.resize(terNum + cellNum);
    ly.resize(terNum + cellNum);
    bw.resize(terNum + cellNum);
    bh.resize(terNum + cellNum);
    idx2name.resize(terNum + cellNum);
    lay.resize(terNum + cellNum);
    olay.resize(terNum + cellNum);
    for (int i=terNum; i<terNum+cellNum; i++)
        fin >> idx2name[i] >> lx[i] >> ly[i] >> bw[i] >> bh[i];

    ox.resize(terNum + cellNum);
    oy.resize(terNum + cellNum);

    fin.close();
}

void output(string& outfile)
{
    ofstream fout(outfile);

    for (int i=terNum; i<terNum+cellNum; i++)
    {
        fout << idx2name[i] << " " 
        << ox[i] << " "
        << oy[i] << " "
        << olay[i] << endl;
    }

    fout.close();
}

void insertTerBlock(int top, int bot, int t)
{
    int l_ = lx[t], r_ = l_ + bw[t];

    for (int i=bot; i<=top; i++)
    {
        bool flag = false;
        for (auto it=ter[i].begin(); it!=ter[i].end(); it++)
        {
            if (r_ < it->first) break;
            if (l_ <= it->second && r_ >= it->first)
            {
                if (l_ < it->first) it->first = l_;
                if (r_ > it->second) it->second = r_;
                flag = true;
            }
        }

        if (!flag) ter[i].emplace_back(make_pair(l_, r_));
    }
}

void updateTerBlocks()
{
    for (int i=0; i<rn; i++)
    {
        for (auto it=ter[i].begin(); it!=ter[i].end(); it++)
        {
            lx0.push_back(it->first);
            lx1.push_back(it->first);
            ly0.push_back(i*rh);
            ly1.push_back(i*rh);
            bw0.push_back(it->second - it->first);
            bw1.push_back(it->second - it->first);
            bh0.push_back(rh);
            bh1.push_back(rh);
        }
        virTerNum += ter[i].size();
    }
}

void buildTwoDieConfig()
{
    lx0.clear();
    lx1.clear();
    ly0.clear();
    ly1.clear();
    bw0.clear();
    bw1.clear();
    bh0.clear();
    bh1.clear();
    d0map.clear();
    d1map.clear();
    virTerNum = 0;
    ter.resize(rn);

    vector<int> ter_order;
    ter_order.resize(terNum);
    for (int i=0; i<terNum; i++)
        ter_order[i] = i;
    sort(ter_order.begin(), ter_order.end(), [&](int a, int b){
        return lx[a] < lx[b];
    });

    for (auto it=ter_order.begin(); it!=ter_order.end(); it++)
    {
        int bot = ly[*it] / rh, top = (ly[*it] + bh[*it]) / rh;

        insertTerBlock(top, bot, *it);
    }
    
    updateTerBlocks();

    for (int i=terNum; i<terNum+cellNum; i++)
    {
        if (lay[i] == 0)
        {
            lx0.push_back(lx[i]);
            ly0.push_back(ly[i]);
            bw0.push_back(bw[i]);
            bh0.push_back(bh[i]);
            d0map.push_back(i);
        }
        else
        {
            lx1.push_back(lx[i]);
            ly1.push_back(ly[i]);
            bw1.push_back(bw[i]);
            bh1.push_back(bh[i]);
            d1map.push_back(i);
        }

        lay[i] = -1;
    }
}

void writeBack(Abacus& die, int l)
{
    int idx; 
    vector<int>& dmap = l == 0 ? d0map : d1map;

    for (int i = virTerNum; i < die.gx.size(); i++)
    {
        idx = dmap[i-virTerNum];
        ox[idx] = die.gx[i];
        oy[idx] = die.gy[i];
        olay[idx] = l;
    }
}

double calCost(Abacus& top, Abacus& bot)
{
    double cost = 0;
    int idx;

    for (int i = virTerNum; i < top.gx.size(); i++)
    {
        idx = d0map[i-virTerNum];
        cost += abs(lx[idx] - top.gx[i]) + abs(ly[idx] - top.gy[i]);
    }

    for (int i = virTerNum; i < bot.gx.size(); i++)
    {
        idx = d1map[i-virTerNum];
        cost += abs(lx[idx] - bot.gx[i]) + abs(ly[idx] - bot.gy[i]);
    }

    return cost;
}

int main(int argc, char* argv[])
{
    timespec bg, tmp, pre;
    int maxT;
    double best_cost = __DBL_MAX__, tmp_cost;
    string in = argv[1], out = argv[2];
    clock_gettime(CLOCK_REALTIME, &bg);

    parse(in);

    do
    {
        clock_gettime(CLOCK_REALTIME, &pre);

        FM fm(terNum, cellNum, W, H, lx, ly, bw, bh, lay);
        fm.start();

        buildTwoDieConfig();

        auto t = std::async(std::launch::async, [&]{
            Abacus top(lx0, ly0, bw0, bh0, W, rh, rn, virTerNum);
            top.start();
            return top;
        });
        
        auto b = std::async(std::launch::async, [&]{
            Abacus bot(lx1, ly1, bw1, bh1, W, rh, rn, virTerNum);
            bot.start();
            return bot;
        });

        auto top = t.get();
        auto bot = b.get();

        tmp_cost = calCost(top, bot);
        if (tmp_cost < best_cost)
        {
            writeBack(top, 0);
            writeBack(bot, 1);
            best_cost = tmp_cost;
        }

        clock_gettime(CLOCK_REALTIME, &tmp);
        if (tmp.tv_sec - pre.tv_sec > maxT) 
            maxT = tmp.tv_sec - pre.tv_sec;

    } while (tmp.tv_sec - bg.tv_sec + maxT < 290);

    output(out);

    return 0;
}