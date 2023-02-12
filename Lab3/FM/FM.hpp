#pragma once

#include <algorithm>
#include <cstdlib> 
#include <vector> 
#include <math.h>
#include <random>
#include <future>
#include <limits.h>
#include <map>
using namespace std;

class FM
{
    public:
        vector<int>& lx;
        vector<int>& ly;
        vector<int>& bw;
        vector<int>& bh;
        vector<int>& lay;
        vector<int> rx;
        vector<int> ry;
        
        vector<int> nds;
        vector<vector<int>> eds;
        vector<vector<int>> wts;
       
        int terNum, cellNum;
        int k, pass;

        void preprocess();
        void initialPlacement(vector<int>&);
        void postprocess(vector<vector<int>>&);
        void start();

        void updateGain(int, int);
        void updateDataStructure(int, vector<int>&, vector<int>&, multimap<int, int>&);
        void revertMoves(vector<int>&, vector<int>&, vector<int>&, multimap<int, int>&);
        void reset(vector<int>&, vector<bool>&, int&, int&);
        int selectVictimNode(vector<int>&, vector<int>&, vector<bool>&, multimap<int, int>&, int&, int&);
        multimap<int, int> buildGainMap(vector<int>&, vector<int>&);
        void Partition(vector<int>&);

        FM(int terNum, int cellNum, int W, int H,
            vector<int>& lx,
            vector<int>& ly,
            vector<int>& bw,
            vector<int>& bh,
            vector<int>& lay): terNum(terNum),
                                cellNum(cellNum),
                                lx(lx),
                                ly(ly),
                                bw(bw),
                                bh(bh),
                                lay(lay),
                                pass(40),
                                k(-1) {
                                    wts.resize(lx.size(), vector<int>({}));
                                    eds.resize(lx.size(), vector<int>({}));
                                };
};

void FM::revertMoves(vector<int>& lay_, vector<int>& revertBuf, vector<int>& gain, multimap<int, int>& gainMap)
{
    for (auto it=revertBuf.rbegin(); it!=revertBuf.rend(); it++)
    {
        updateDataStructure(*it, lay_, gain, gainMap);
    }
}

void FM::updateDataStructure(int n_idx, vector<int>& lay_, vector<int>& gain, multimap<int, int>& gainMap)
{
    int delta, to;
    vector<int> pre_g = {gain[n_idx]};

    for (int i=0; i<eds[n_idx].size(); i++)
    {
        delta = 2*wts[n_idx][i];
        to = eds[n_idx][i];
        pre_g.push_back(gain[to]);

        if (lay_[n_idx] == lay_[to])
        {
            gain[n_idx] -= delta;
            gain[to] -= delta;
        }
        else
        {
            gain[n_idx] += delta;
            gain[to] += delta;
        }
    }

    auto itr = gainMap.find(pre_g[0]);
    for (int i=0; i<gainMap.count(pre_g[0]); i++, itr++)
    {
        if (itr->second == n_idx)
        {
            gainMap.erase(itr);
            gainMap.insert({gain[n_idx], n_idx});
            break;
        }
    }

    for (int i=0; i<eds[n_idx].size(); i++)
    {
        int idx = eds[n_idx][i];
        itr = gainMap.find(pre_g[i+1]);

        for (int j=0; j<gainMap.count(pre_g[i+1]); j++, itr++)
        {
            if (itr->second == idx)
            {
                gainMap.erase(itr);
                gainMap.insert({gain[idx], idx});
                break;
            }
        }
    }

    lay_[n_idx] = 1 - lay_[n_idx];
}

int FM::selectVictimNode(vector<int>& lay_,
                         vector<int>& revertBuf,
                         vector<bool>& locked,
                         multimap<int, int>& gainMap,
                         int& debt,
                         int& noImprove)
{
    int n_idx = -1, g;

    for (auto it=gainMap.rbegin(); it!=gainMap.rend(); it++)
    {
        if (!locked[it->second])
        {
            locked[it->second] = true;
            n_idx = it->second;
            g = it->first;
            break;
        }
    }

    if (n_idx != -1)
    {
        debt += g;

        if (g < 0) noImprove++;

        if (debt < 0)
        {
            revertBuf.push_back(n_idx);
            // noImprove++;
        }
        else
        {
            revertBuf.clear();
            // noImprove = 0;
            debt = 0;
        }

        if (noImprove == 100) {n_idx = -1;}
        // if (noImprove > k) {n_idx = -1;}
    }

    return n_idx;
}

void FM::reset(vector<int>& revertBuf, vector<bool>& locked, int& noImprove, int& debt)
{
    revertBuf.clear();
    noImprove = 0;
    debt = 0;

    for (int i=0; i<locked.size(); i++)
        locked[i] = false;
}

multimap<int, int> FM::buildGainMap(vector<int>& lay_, vector<int>& gain)
{
    multimap<int, int> gainMap;

    for (auto it=nds.begin(); it!=nds.end(); it++)
    {
        auto &neigs = eds[*it];

        for (int i=0; i<neigs.size(); i++)
        {
            if (lay_[*it] == lay_[neigs[i]])
                gain[*it] += wts[*it][i];
            else
                gain[*it] -= wts[*it][i];
        }
    }

    for (auto it=nds.begin(); it!=nds.end(); it++)
        gainMap.insert({gain[*it], *it});

    return gainMap;
}

void FM::Partition(vector<int>& lay_){
    vector<int> gain(lx.size(), 0), revertBuf;
    auto gainMap = buildGainMap(lay_, gain);
    vector<bool> locked(lx.size(), false);
    int noImprove = 0, debt = 0;

    // implement FM algo
    for (int i=0; i<pass; i++)
    {
        while (true)
        {
            int n_idx = selectVictimNode(lay_, revertBuf, locked, gainMap, debt, noImprove);

            // no available node found
            if (n_idx == -1) break;

            // update gain to neighboring nodes
            updateDataStructure(n_idx, lay_, gain, gainMap);
        }

        // revert previous moves to the smallest cut size
        if (debt < 0) 
        {
            revertMoves(lay_, revertBuf, gain, gainMap);
        }

        reset(revertBuf, locked, noImprove, debt);
    }
}


void FM::preprocess()
{
    rx.resize(terNum+cellNum);
    ry.resize(terNum+cellNum);
    for (int i=terNum; i<terNum+cellNum; i++)
    {
        rx[i] = lx[i] + bw[i];
        ry[i] = ly[i] + bh[i];
    }

    vector<int> order;
    order.reserve(cellNum);
    for (int i=terNum; i<terNum+cellNum; i++)
        order.push_back(i);
    sort(order.begin(), order.end(), [&](int a, int b){
        return lx[a] < lx[b];
    });

    vector<bool> picked(lx.size(), false);
    int a, b, weight = -1;
    for (int i=0; i<cellNum; i++)
    {
        bool flag = false;
        a = order[i];
        for (int j=i+1; j<cellNum; j++)
        {
            b = order[j];
            if (lx[a] < rx[b] && lx[b] < rx[a] && ly[b] < ry[a] && ly[a] < ry[b])
            {
                weight = (min(rx[a], rx[b]) - max(lx[a], lx[b])) *
                            (min(ry[a], ry[b]) - max(ly[a], ly[b]));
                wts[a].push_back(weight);
                wts[b].push_back(weight);
                eds[a].push_back(b);
                eds[b].push_back(a);
                flag = true;

                if (!picked[b]) 
                {
                    nds.push_back(b);
                    picked[b] = true;
                }
            }

            if (rx[a] <= lx[b]) break;
        }

        if (flag && !picked[a]) 
        {
            nds.push_back(a);
            picked[a] = true;
        }
    }

    k = nds.size() * 0.1;
}

void FM::initialPlacement(vector<int>& lay_)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<bool> nodesPicked(lx.size(), false);
    std::vector<int> order;
    order.reserve(nds.size());
    for (int i=0; i<nds.size(); i++)
        order.push_back(nds[i]);
    std::shuffle(order.begin(), order.end(), gen);

    for (auto it=order.begin(); it!=order.end(); it++)
        lay_[*it] = gen()%2;
}

void FM::postprocess(vector<vector<int>>& lays)
{
    int best_score = 0, idx, a;

    for (int i=0; i<lays.size(); i++)
    {
        int tmp_score = 0;
        auto& lay_ = lays[i];

        for (int j=0; j<nds.size(); j++)
        {
            a = nds[j];
            for (int k=0; k<eds[a].size(); k++)
            {
                if (lay_[a] != lay_[eds[a][k]])
                    tmp_score += wts[a][k];
            }
        }

        if (tmp_score > best_score)
        {
            idx = i;
            best_score = tmp_score;
        }
    }

    auto& lay_ = lays[idx];
    int p0size = 0, p1size = 0;
    for (auto it=nds.begin(); it!=nds.end(); it++)
    {
        lay[*it] = lay_[*it];

        if (lay_[*it] == 0)
            p0size += (rx[*it] - lx[*it]) * (ry[*it] - ly[*it]);
        else
            p1size += (rx[*it] - lx[*it]) * (ry[*it] - ly[*it]);
    }

    for (int i=terNum; i<terNum+cellNum; i++)
    {
        if (lay[i] != -1) continue;
        if (p0size > p1size)
        {
            lay[i] = 1;
            p1size += (rx[i] - lx[i]) * (ry[i] - ly[i]);
        }
        else
        {
            lay[i] = 0;
            p0size += (rx[i] - lx[i]) * (ry[i] - ly[i]);
        }
    }
}

void FM::start()
{
    const int round = 10;
    vector<vector<int>> lays(round, vector<int>(lay.size(), -1));
    vector<future<void>> futs;

    preprocess();

    for (int i=0; i<round; i++)
    {        
        futs.push_back(std::async(std::launch::async, [&](vector<int>& lay_){
            initialPlacement(lay_);
            Partition(lay_);
        }, std::ref(lays[i])));
    }

    for (int i=0; i<round; i++)
        futs[i].wait();

    postprocess(lays);
}