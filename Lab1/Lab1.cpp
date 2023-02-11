#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <memory>
#include <cstring>
#include <assert.h>
#include <algorithm>

using namespace std;

#define nullPtr shared_ptr<block>(nullptr)
#define draw false

struct block
{
    shared_ptr<block> rt, tr, bl, lb;
    int idx, x0, y0, x1, y1;

    block(int x0, int y0, int x1, int y1): idx(-1), x0(x0), y0(y0), x1(x1), y1(y1),
                                    rt(nullPtr), tr(nullPtr), bl(nullPtr), lb(nullPtr) {};
};


// split the given block into 2 blocks vertically
// original block on the bottom, new generated block on top
// y is the y-coordinate of the new generated block
// return the new generated block (top)
shared_ptr<block> splitY(shared_ptr<block> tmp, int y)
{
    assert(tmp->y0 < y && tmp->y1 > y);

    auto newBlk = make_shared<block>(tmp->x0, y, tmp->x1, tmp->y1);

    // top adjaceny
    for (auto blk = tmp->rt; blk != nullPtr && blk->x0 >= tmp->x0; blk = blk->bl)
        blk->lb = newBlk;
    newBlk->rt = tmp->rt;

    // bottom adjaceny
    tmp->rt = newBlk;
    newBlk->lb = tmp;

    // right adjaceny
    auto blk = tmp->tr;
    for (; blk != nullPtr && blk->y0 >= y; blk = blk->lb)
        blk->bl = newBlk;
    newBlk->tr = tmp->tr;
    tmp->tr = blk;
    
    // left adjaceny
    blk = tmp->bl;
    while (blk != nullPtr && blk->y1 <= y)
        blk = blk->rt;
    newBlk->bl = blk;
    for (; blk != nullPtr && blk->y1 <= tmp->y1; blk = blk->rt)
        blk->tr = newBlk;

    // resize original block
    tmp->y1 = y;

    return newBlk;
}

// split the given into 3 blocks horizontally
// orginal block on the left, midBlock in the middle, rightBlk ont the right
// x1 is the x0 value of midBlk, x2 is the x0 value of rightBlk
// return the middle block
shared_ptr<block> split3X(shared_ptr<block> tmp, int x1, int x2)
{
    int right_boundary = tmp->x1, left_boundary = tmp->x0;
    auto leftBlk = tmp->x0 == x1 ? tmp : make_shared<block>(left_boundary, tmp->y0, x1, tmp->y1);
    auto rightBlk = tmp->x1 == x2 ? tmp : make_shared<block>(x2, tmp->y0, right_boundary, tmp->y1);

    // top adjaceny
    auto blk = tmp->rt;
    rightBlk->rt = blk;
    for (; blk != nullPtr && blk->x0 >= x2; blk = blk->bl)
        blk->lb = rightBlk;   
    tmp->rt = blk;
    for (; blk != nullPtr && blk->x0 >= x1; blk = blk->bl) {}
    if (tmp != leftBlk) 
    {
        leftBlk->rt = blk;
        for (; blk != nullPtr && blk->x0 >= left_boundary; blk = blk->bl)
            blk->lb = leftBlk;
    }

    // bottom adjaceny
    blk = tmp->lb;
    leftBlk->lb = blk;
    for (; blk != nullPtr && blk->x1 <= x1; blk = blk->tr) 
        blk->rt = leftBlk;
    tmp->lb = blk;
    for (; blk != nullPtr && blk->x1 <= x2; blk = blk->tr) {}
    if (tmp != rightBlk) 
    {
        rightBlk->lb = blk;
        for (; blk != nullPtr && blk->x1 <= right_boundary; blk = blk->tr)
            blk->rt = rightBlk; 
    }

    // left and right adjacency
    if (leftBlk != tmp)
    {
        // left adjaceny
        for (auto blk = tmp->bl; blk != nullPtr && blk->y1 <= tmp->y1; blk = blk->rt)
            blk->tr = leftBlk;

        leftBlk->tr = tmp;
        leftBlk->bl = tmp->bl;
        tmp->bl = leftBlk;
        tmp->x0 = x1;
    }
    if (rightBlk != tmp)
    {
        // right adjaceny
        for (auto blk = tmp->tr; blk != nullPtr && blk->y0 >= tmp->y0; blk = blk->lb)
            blk->bl = rightBlk;

        rightBlk->bl = tmp;
        rightBlk->tr = tmp->tr;
        tmp->tr = rightBlk;
        tmp->x1 = x2;
    }

    return tmp;
}

// merge block b2 to b1 vertically
// return the merged block
shared_ptr<block> mergeY(shared_ptr<block> b1, shared_ptr<block> b2)
{
    // swap b1 and b2 if b2 is lower than b1
    if (b2->y0 < b1->y0)
    {
        auto tmp = b1;
        b1 = b2;
        b2= tmp;
    }

    // top adjacency
    for (auto blk = b2->rt; blk != nullPtr && blk->x0 >= b2->x0; blk = blk->bl)
        blk->lb = b1;
    b1->rt = b2->rt;

    // right adjacency
    for (auto blk = b2->tr; blk != nullPtr && blk->y0 >= b2->y0; blk = blk->lb)
        blk->bl = b1;
    b1->tr = b2->tr;

    // left adjacency
    for (auto blk = b2->bl; blk != nullPtr && blk->y1 <= b2->y1; blk = blk->rt)
        blk->tr = b1;

    // resize b1
    b1->y1 = b2->y1;

    return b1;
}

// return block containing the point with coordinate (x, y)
shared_ptr<block> find_point(shared_ptr<block> tmp, int x, int y)
{
    while (tmp->x0 > x || tmp->x1 <= x || tmp->y0 > y || tmp->y1 <= y)
    {
        assert(tmp != nullPtr);

        while (tmp->y1 <= y)
            tmp = tmp->rt;

        while (tmp->y0 > y)
            tmp = tmp->lb;

        while (tmp->x0 > x)
            tmp = tmp->bl;

        while (tmp->x1 <= x)
            tmp = tmp->tr;
    }

    return tmp;
}

void block_insert(shared_ptr<block> base, vector<string> info)
{
    int bottom = stoi(info[2]), 
        top = bottom + stoi(info[4]), 
        left = stoi(info[1]), 
        right = left + stoi(info[3]);

    // split block horizontally on the top boundary
    // auto top_block = find_point(base, right, top);
    auto top_block = find_point(base, left, top-1);
    if (top_block->y1 != top)
    {
        top_block = splitY(top_block, top);
        if (top_block->rt != nullPtr && 
            top_block->rt->idx == -1 &&
            top_block->rt->x0 == top_block->x0 && 
            top_block->rt->x1 == top_block->x1)
            mergeY(top_block->rt, top_block);
    }

    // split block horizontally on the lower boundary
    auto tmp = find_point(base, left, bottom);
    if (tmp->y0 != bottom)
    {
        tmp = splitY(tmp, bottom);
        auto bottom_nei = tmp->lb;
        if (bottom_nei->lb != nullPtr && 
            bottom_nei->lb->idx == -1 &&
            bottom_nei->lb->x0 == bottom_nei->x0 && 
            bottom_nei->lb->x1 == bottom_nei->x1)
            mergeY(bottom_nei->lb, bottom_nei);
    }

    // split blocks into three pieces vertically from top to bottom
    // merge blocks inside target region vertically
    // merge left/right neighbors if needed
    while (true)
    {
        tmp = split3X(tmp, left, right);

        if (tmp->y0 != bottom)
        {
            auto right_nei = tmp->tr;
            if (right_nei != nullPtr && 
                right_nei->lb != nullPtr &&
                right_nei->x1 == right_nei->lb->x1 && 
                right_nei->idx == right_nei->lb->idx)
                mergeY(right_nei->lb, right_nei);

            auto left_nei = tmp->bl;
            if (left_nei != nullPtr && 
                left_nei->lb != nullPtr &&
                left_nei->x0 == left_nei->lb->x0 &&
                left_nei->x1 == left_nei->lb->x1 && 
                left_nei->idx == left_nei->lb->idx)
                mergeY(left_nei->lb, left_nei);

            tmp = mergeY(tmp->lb, tmp);
        }

        if (tmp->y1 == top)
        {
            tmp->idx = stoi(info[0]);
            break;
        }

        tmp = tmp->rt;
    }
}

// enumerate for drawing
void enumerate_d(shared_ptr<block> tmp, vector<vector<int>>& info, int& cnt, int bottom)
{
    while (tmp != nullPtr && tmp->y1 > bottom)
    {   
        info.emplace_back(vector<int>({tmp->idx, 
                                        tmp->x0, 
                                        tmp->y0, 
                                        (tmp->x1 - tmp->x0), 
                                        (tmp->y1 - tmp->y0)}));

        if (tmp->tr != nullPtr && tmp->tr->y0 >= tmp->y0)
            enumerate_d(tmp->tr, info, cnt, tmp->y0);

        tmp = tmp->lb;
        cnt++;
    }
}

// return number of neighboring block and space
vector<int> find_block_and_space(shared_ptr<block> tmp)
{
    vector<int> ans({tmp->idx, 0, 0});

    // top adjacency
    for (auto blk = tmp->rt; blk != nullPtr && blk->x1 > tmp->x0; blk = blk->bl)
        blk->idx == -1 ? ans[2]++ : ans[1]++;

    // bottom adjacency
    for (auto blk = tmp->lb; blk != nullPtr && blk->x0 < tmp->x1; blk = blk->tr)
        blk->idx == -1 ? ans[2]++ : ans[1]++;

    // left adjacency
    for (auto blk = tmp->bl; blk != nullPtr && blk->y0 < tmp->y1; blk = blk->rt)
        blk->idx == -1 ? ans[2]++ : ans[1]++;

    // right adjacency
    for (auto blk = tmp->tr; blk != nullPtr && blk->y1 > tmp->y0; blk = blk->lb)
        blk->idx == -1 ? ans[2]++ : ans[1]++;

    return ans;
}

// enumerate for output
void enumerate_o(shared_ptr<block> tmp, vector<vector<int>>& info, int& cnt, int bottom)
{
    while (tmp != nullPtr && tmp->y0 >= bottom)
    {   
        if (tmp->idx != -1)
            info.emplace_back(find_block_and_space(tmp));

        if (tmp->tr != nullPtr && tmp->tr->y0 >= tmp->y0)
            enumerate_o(tmp->tr, info, cnt, tmp->y0);

        tmp = tmp->lb;
        cnt++;
    }
}

// write to layout.txt for drawing
void enumerate_and_output(shared_ptr<block> base, int x, int y, string out)
{
    auto tmp = find_point(base, 0, y-1);
    vector<vector<int>> info;
    ofstream fout;
    int cnt = 0;

    if (draw)
    {
        enumerate_d(tmp, info, cnt, 0);

        fout.open("layout.txt");
        fout << cnt << endl;
        fout << x << " " << y << endl;

        for (auto const &it : info)
        fout << it[0] << " " << it[1] << " " << it[2] << " " << it[3] << " " << it[4] << endl;
    }
    else
    {
        enumerate_o(tmp, info, cnt, 0);

        sort(info.begin(), info.end(), [&](vector<int>& a, vector<int>& b) {
            return a[0] < b[0];
        });

        fout.open(out);
        fout << cnt << endl;

        for (auto const &it : info)
            fout << it[0] << " " << it[1] << " " << it[2] << endl;
    }

    fout.close();
}

int main(int argc, char *argv[])
{
    ifstream fin(argv[1]);

    string line;
    int width, height;

    fin >> width >> height;
    fin.ignore(numeric_limits<std::streamsize>::max(), '\n');

    auto base = make_shared<block>(0, 0, width, height);
    vector<vector<int>> find_buf;

    while(getline(fin, line))
    {
        vector<string> lineVec;
        stringstream s_str(line);
        string str;

        while(getline(s_str, str, ' '))
        {
            lineVec.push_back(str);
        }

        if (lineVec[0] == "P")
        {
            auto blk = find_point(base, stoi(lineVec[1]), stoi(lineVec[2]));
            find_buf.emplace_back(vector<int>({blk->x0, blk->y0}));
        }
        else
        {
            block_insert(base, lineVec);
        }
    }

    fin.close();
    enumerate_and_output(base, width, height, argv[2]);
    
    if (!draw)
    {
        ofstream fout(argv[2], std::ios_base::app);

        for (auto const &it : find_buf)
            fout << it[0] << " " << it[1] << endl;

        fout.close();
    }

    return 0;
}