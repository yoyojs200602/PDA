#include <vector>
#include <algorithm>
#include <limits.h>

using namespace std;


class Row
{
    public:
        int y, xmin, xmax, totalW;
        vector<int> cell;

        Row(int y, int X, int W): y(y), xmin(X), xmax(W), totalW(0), cell({}) {};
};

class Cluster
{
    public:
        int x, w, q, e;
        int nFirst, nLast;

        Cluster(): x(-1), w(0), q(0), e(0), nFirst(-1), nLast(-1) {};
};

class Abacus
{
    public:
        vector<Row> r;
        vector<int>& gx;
        vector<int>& gy;
        vector<int>& w;
        vector<int>& h;
        vector<vector<int>> ter;

        int W, rh, rn, terNum;

        Abacus(vector<int>& x,
                vector<int>& y,
                vector<int>& w,
                vector<int>& h,
                int W,
                int rh,
                int rn_,
                int terNum
            ): W(W), gx(x), gy(y), w(w), h(h), rh(rh), terNum(terNum) {
                rn = 0;
                ter.resize(rn_);
                for (int i=0; i<terNum; i++) {ter[y[i]/rh].push_back(i);}

                for (int i=0; i<rn_; i++) {
                    int l_ = 0, h_ = rh*i;
                    
                    for (auto it=ter[i].begin(); it!=ter[i].end(); it++)
                    {
                        if (x[*it] != 0)
                        {
                            r.emplace_back(Row(h_, l_, x[*it]));
                            rn++;
                        }

                        l_ = x[*it] + w[*it];
                    }

                    if (l_ != W)
                    {
                        r.emplace_back(Row(h_, l_, W));
                        rn++;
                    }
                }
            };

        bool start();
        int findRow(int);
        int placeRow(int, bool);
        void addCell(Cluster&, int);
        void addCluster(Cluster&, Cluster&);
        void collapse(vector<Cluster>&, int);
};


int Abacus::findRow(int idx)
{
    int top = rn-1, bot = 0, mid, idx_y = gy[idx];

    while (top - bot > 1)
    {
        mid = (top + bot) >> 1;
        r[mid].y > idx_y ? top = mid : bot = mid;
    }

    return bot;
}

int Abacus::placeRow(int r_idx, bool flag)
{
    auto& cells = r[r_idx].cell;

    if (r[r_idx].totalW + w[cells.back()] > r[r_idx].xmax - r[r_idx].xmin) return INT_MAX;

    vector<Cluster> clv;

    for (auto it=cells.begin(); it!=cells.end(); it++)
    {
        int x_ = gx[*it];
        if  (x_ < r[r_idx].xmin) gx[*it] = r[r_idx].xmin;
        if  (x_ > r[r_idx].xmax - w[*it]) gx[*it] = r[r_idx].xmax - w[*it];
        
        if (clv.size() == 0 || gx[*it] > clv.back().x + clv.back().w)
        {
            clv.emplace_back(Cluster());
            clv.back().x = gx[*it];
            clv.back().nFirst = *it;
            addCell(clv.back(), *it);
        }
        else
        {
            addCell(clv.back(), *it);
            collapse(clv, r_idx);
        }

        gx[*it] = x_;
    }

    if (flag)
    {
        int i = 0, x_ = r[r_idx].xmin;
        for (auto it=clv.begin(); it!=clv.end(); it++)
        {
            x_ = it->x;

            while (true)
            {
                gx[cells[i]] = x_;
                x_ += w[cells[i]];
                if (cells[i++] == it->nLast) break;
            }
        }
    }

    return abs(gx[cells.back()] - (clv.back().x + clv.back().w - w[cells.back()])) + abs(gy[cells.back()] - r[r_idx].y);
}

void Abacus::addCell(Cluster& c, int i)
{
    c.nLast = i;
    c.e++;
    c.q += (gx[i] - c.w);
    c.w += w[i];
}

void Abacus::addCluster(Cluster& c1, Cluster& c2)
{
    c1.nLast = c2.nLast;
    c1.e += c2.e;
    c1.q += (c2.q - c2.e * c1.w);
    c1.w += c2.w;
}

void Abacus::collapse(vector<Cluster>& clv, int r_idx)
{
    int i = clv.size()-1;

    while (true)
    {
        auto& tmp = clv[i];
        tmp.x = tmp.q / tmp.e;

        if (tmp.x < r[r_idx].xmin) {tmp.x = r[r_idx].xmin;}
        if (tmp.x > r[r_idx].xmax - tmp.w) {tmp.x = r[r_idx].xmax - tmp.w;}

        if (i != 0 && clv[i-1].x + clv[i-1].w > tmp.x)
        {
            addCluster(clv[i-1], tmp);
            clv.pop_back();
        }
        else
        {
            break;
        }

        i--;
    }
}

bool Abacus::start()
{
    vector<int> x_order;
    x_order.reserve(gx.size() - terNum);
    for (int i=terNum; i<gx.size(); i++) {x_order.push_back(i);}
    sort(x_order.begin(), x_order.end(), [&](int a, int b){
        return gx[a] < gx[b];
    });

    for (auto it=x_order.begin(); it!=x_order.end(); it++)
    {
        int best_cost = INT_MAX, tmp_cost;
        int rbest = -1, origin = findRow(*it);

        for (int up = origin+1; up < rn; up++)
        {
            if (r[up].y - gy[*it] > best_cost) break;

            r[up].cell.push_back(*it);
            tmp_cost = placeRow(up, false);
            
            if (tmp_cost < best_cost)
            {
                best_cost = tmp_cost;
                rbest = up;
            }

            r[up].cell.pop_back();
        }

        for (int bot = origin; bot >= 0; bot--)
        {
            if (gy[*it] - r[bot].y > best_cost) break;

            r[bot].cell.push_back(*it);
            tmp_cost = placeRow(bot, false);
            
            if (tmp_cost < best_cost)
            {
                best_cost = tmp_cost;
                rbest = bot;
            }

            r[bot].cell.pop_back();
        }

        if (rbest == -1)
            return false;

        if (gx[*it] < r[rbest].xmin) gx[*it] = r[rbest].xmin;
        if (gx[*it] > r[rbest].xmax - w[*it]) gx[*it] = r[rbest].xmax - w[*it];
        r[rbest].cell.push_back(*it);
        placeRow(rbest, true);
        gy[*it] = r[rbest].y;
        r[rbest].totalW += w[*it];
    }

    return true;
}