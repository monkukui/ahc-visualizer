#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("Ofast")

#include "atcoder/modint.hpp"
#include "atcoder/dsu.hpp"
#include <random>
#include <iostream>
#include <cstdio>
#include <string>
#include <regex>
#include <cstring>
#include <deque>
#include <list>
#include <queue>
#include <unordered_set>
#include <stack>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include <set>
#include <complex>
#include <cmath>
#include <limits>
#include <climits>
#include <ctime>
#include <cassert>
#include <numeric>
#include <functional>
#include <bitset>
#include <cstddef>
#include <type_traits>
#include <vector>
#include <sys/time.h>

using namespace std;

using lint = long long int;
using i128 = __int128;
const long long int INF = numeric_limits<long long int>::max() / 4;
const int inf = 1e9;
const long long int MOD = 998244353;
const double MATH_PI = 3.1415926535897932;

bool out_of_grid(int i, int j, int h, int w) {
    if (i < 0) return true;
    if (j < 0) return true;
    if (i >= h) return true;
    if (j >= w) return true;
    return false;
}

template<typename T1, typename T2>
inline void chmin(T1 &a, const T2 &b) { if (a > b) a = b; }

template<typename T1, typename T2>
inline void chmax(T1 &a, const T2 &b) { if (a < b) a = b; }

#define ALL(a) a.begin(),a.end()
#define RALL(a) a.rbegin(),a.rend()
#define rep(i, n) for(int i=0;i<(int)(n);i++)

#define SUM(v) accumulate(ALL(v), 0LL)
#define MIN(v) *min_element(ALL(v))
#define MAX(v) *max_element(ALL(v))

#ifdef LOCAL

ostream& operator<<(ostream& os, __int128_t x) {
    if (x < 0) {
        os << "-";
        x *= -1;
    }
    if (x == 0) {
        return os << "0";
    }
    string s;
    while (x) {
        s += char(x % 10 + '0');
        x /= 10;
    }
    reverse(s.begin(), s.end());
    return os << s;
}
ostream& operator<<(ostream& os, __uint128_t x) {
    if (x == 0) {
        return os << "0";
    }
    string s;
    while (x) {
        s += char(x % 10 + '0');
        x /= 10;
    }
    reverse(s.begin(), s.end());
    return os << s;
}

template <class T, class U>
ostream& operator<<(ostream& os, const pair<T, U>& p);
template <class T> ostream& operator<<(ostream& os, const vector<T>& v);
template <class T> ostream& operator<<(ostream& os, const deque<T>& v);
template <class T, size_t N>
ostream& operator<<(ostream& os, const array<T, N>& a);
template <class T> ostream& operator<<(ostream& os, const set<T>& s);
template <class T, class U>
ostream& operator<<(ostream& os, const map<T, U>& m);

template <typename T, typename Container>
ostream &operator<<(ostream &os,
                         const priority_queue<T, Container> &pq);

template <class T, class U>
ostream& operator<<(ostream& os, const pair<T, U>& p) {
    return os << "P(" << p.first << ", " << p.second << ")";
}

template <class T> ostream& operator<<(ostream& os, const vector<T>& v) {
    os << "[";
    bool f = false;
    for (auto d : v) {
        if (f) os << ", ";
        f = true;
        os << d;
    }
    return os << "]";
}

template <class T> ostream& operator<<(ostream& os, const deque<T>& v) {
    os << "[";
    bool f = false;
    for (auto d : v) {
        if (f) os << ", ";
        f = true;
        os << d;
    }
    return os << "]";
}
template <class T, size_t N>
ostream& operator<<(ostream& os, const array<T, N>& a) {
    os << "[";
    bool f = false;
    for (auto d : a) {
        if (f) os << ", ";
        f = true;
        os << d;
    }
    return os << "]";
}

template <class T> ostream& operator<<(ostream& os, const set<T>& s) {
    os << "{";
    bool f = false;
    for (auto d : s) {
        if (f) os << ", ";
        f = true;
        os << d;
    }
    return os << "}";
}
template <class T> ostream& operator<<(ostream& os, const multiset<T>& s) {
    os << "{";
    bool f = false;
    for (auto d : s) {
        if (f) os << ", ";
        f = true;
        os << d;
    }
    return os << "}";
}

template <class T, class U>
ostream& operator<<(ostream& os, const map<T, U>& s) {
    os << "{";
    bool f = false;
    for (auto p : s) {
        if (f) os << ", ";
        f = true;
        os << p.first << ": " << p.second;
    }
    return os << "}";
}

template <typename T, typename Container>
ostream &operator<<(ostream &os,
                         const priority_queue<T, Container> &pq) {
  priority_queue<T, Container> temp =
      pq;

  os << "[";
  bool f = false;
  while (!temp.empty()) {
    if (f) os << ", ";
    f = true;
    os << temp.top();
    temp.pop();
  }
  os << "]";

  return os;
}

struct PrettyOS {
    ostream& os;
    bool first;

    template <class T> auto operator<<(T&& x) {
        if (!first) os << ", ";
        first = false;
        os << x;
        return *this;
    }
};
template <class... T> void dbg0(T&&... t) {
    (PrettyOS{cerr, true} << ... << t);
}
#define dbg(...)                                            \
    do {                                                    \
        cerr << __LINE__ << " : " << #__VA_ARGS__ << " = "; \
        dbg0(__VA_ARGS__);                                  \
        cerr << endl;                                       \
    } while (false);
#else
#define dbg(...)
#endif

namespace monkukui {
    namespace output {
        void yesno(bool f) {
            if (f) cout << "Yes" << endl;
            else cout << "No" << endl;
        }

        template<class T>
        void vec(vector <T> a) {
            for (int i = 0; i < (int) a.size(); i++) {
                if (i + 1 == (int) a.size()) cout << a[i] << endl;
                else cout << a[i] << " ";
            }
        }

        template<class T>
        void grid(vector <vector<T>> g) {
            int h = g.size();
            int w = g[0].size();
            for (int i = 0; i < h; i++) {
                for (int j = 0; j < w; j++) {
                    if (j + 1 == w) cout << g[i][j] << endl;
                    else cout << g[i][j] << " ";
                }
            }
        }
    }

    namespace input {
        template<class T>
        vector <T> vec(int n) {
            vector <T> a(n);
            for (int i = 0; i < n; i++) cin >> a[i];
            return a;
        }

        vector <vector<int>> unweighted_graph(int n, int m, bool directed = false) {
            vector <vector<int>> g(n);
            for (int i = 0; i < m; i++) {
                int a, b;
                cin >> a >> b;
                a--;
                b--;
                g[a].emplace_back(b);
                if (!directed) g[b].emplace_back(a);
            }
            return g;
        }

        template<class T>
        vector <vector<pair < int, T>>>

        weighted_graph(int n, int m, bool directed = false) {
            vector <vector<int>> g(n);
            for (int i = 0; i < m; i++) {
                int a, b;
                cin >> a >> b;
                T c;
                cin >> c;
                a--;
                b--;
                g[a].emplace_back({b, c});
                if (!directed) g[b].emplace_back({a, c});
            }
            return g;
        }

        template<class T>
        vector <vector<T>> grid(int h, int w) {
            vector <vector<T>> ret(h, vector<T>(w));
            for (int i = 0; i < h; i++) {
                for (int j = 0; j < w; j++) {
                    cin >> ret[i][j];
                }
            }
            return ret;
        }
    }
}

using mint = atcoder::modint998244353;

struct Timer {
    struct timeval start, cur;
    double limit;

    Timer() : limit(0) { gettimeofday(&start, NULL); }

    Timer(double l) : limit(l) { gettimeofday(&start, NULL); }

    bool isLimit() { return curTime() > limit; }

    double curTime() {
        gettimeofday(&cur, NULL);
        return (cur.tv_sec - start.tv_sec) + (cur.tv_usec - start.tv_usec) / 1e6;
    }
};

// xrand.h - An implementation of xorshift random number generator.
// See the description of class XRand for its usage.
// Version 0.2

// Copyright (c) 2013  Kazuhiro Hosaka (hos@hos.ac)
//
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the authors be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
//    1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
//
//    2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
//
//    3. This notice may not be removed or altered from any source
//    distribution.

#ifndef XRAND_H_
#define XRAND_H_

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <algorithm>

#ifndef UINT32_MAX
#define UINT32_MAX (4294967295U)
#endif

#ifndef UINT64_MAX
#define UINT64_MAX (18446744073709551615ULL)
#endif

#define rand           _DO_NOT_USE_RAND_
#define srand          _DO_NOT_USE_SRAND_
#define random_shuffle _DO_NOT_USE_RANDOM_SHUFFLE_

class XRand {
public:
    // Initializes with the given seed.
    explicit XRand(uint64_t seed = 0);

    void Reset(uint64_t seed);

    // Generates the next random number, changing the state of the generator.
    // NextULong() calls NextUInt() twice.
    uint32_t NextUInt();

    uint64_t NextULong();

    // Returns an integer, almost uniformly distributed in [0, m).
    // It must hold that m != 0.
    // Calls NextUInt() once for UInt and twice for ULong.
    uint32_t NextUInt(uint32_t m);

    uint64_t NextULong(uint64_t m);

    // Returns an integer, almost uniformly distributed in [a, b].
    // It must hold that a <= b.
    // Calls NextUInt() once for Int and twice for Long.
    int32_t NextInt(int32_t a, int32_t b);

    int64_t NextLong(int64_t a, int64_t b);

    // The unbiased version of the methods above.
    // The number of calls to NextUInt() is not deterministic.
    // The expected number of that is at most two for (U)Int and four for (U)Long.
    uint32_t NextUIntUnbiased(uint32_t m);

    uint64_t NextULongUnbiased(uint64_t m);

    int32_t NextIntUnbiased(int32_t a, int32_t b);

    int64_t NextLongUnbiased(int64_t a, int64_t b);

    // Returns a double, uniformly distributed in [0.0, 1.0).
    // Calls NextUInt() twice.
    double NextDouble();

    // Returns a double, normally distributed with expectation 0 and variance 1.
    // Calls NextUInt() four times.
    double NextGaussian();

    // Shuffles the interval [first, last).
    // Calls NextUInt() (last - first) times.
    template<typename RandomAccessIterator>
    void Shuffle(RandomAccessIterator first, RandomAccessIterator last);

private:
    uint32_t x_, y_, z_, w_;
};

XRand::XRand(uint64_t seed) {
    Reset(seed);
}

void XRand::Reset(uint64_t seed) {
    x_ = 314159265;
    y_ = 358979323;
    z_ = 846264338 ^ (65535 * static_cast<uint32_t>(seed >> 32));
    w_ = 327950288 ^ (65535 * static_cast<uint32_t>(seed));
}

uint32_t XRand::NextUInt() {
    const uint32_t t = x_ ^ x_ << 11;
    x_ = y_;
    y_ = z_;
    z_ = w_;
    return w_ = w_ ^ w_ >> 19 ^ t ^ t >> 8;
}

uint64_t XRand::NextULong() {
    const uint64_t high = NextUInt();
    const uint64_t low = NextUInt();
    return high << 32 | low;
}

uint32_t XRand::NextUInt(uint32_t m) {
    assert(m != 0);
    return NextUInt() % m;
}

uint64_t XRand::NextULong(uint64_t m) {
    assert(m != 0);
    return NextULong() % m;
}

int32_t XRand::NextInt(int32_t a, int32_t b) {
    assert(a <= b);
    return a + NextUInt(b - a + 1);
}

int64_t XRand::NextLong(int64_t a, int64_t b) {
    assert(a <= b);
    return a + NextULong(b - a + 1);
}

uint32_t XRand::NextUIntUnbiased(uint32_t m) {
    assert(m != 0);
    if (m & (m - 1)) {
        const uint32_t limit = UINT32_MAX / m * m;
        for (;;) {
            const uint32_t value = NextUInt();
            if (value < limit) {
                return value % m;
            }
        }
    } else {
        return NextUInt() & (m - 1);
    }
}

uint64_t XRand::NextULongUnbiased(uint64_t m) {
    assert(m != 0);
    if (m & (m - 1)) {
        const uint64_t limit = UINT64_MAX / m * m;
        for (;;) {
            const uint64_t value = NextULong();
            if (value < limit) {
                return value % m;
            }
        }
    } else {
        return NextULong() & (m - 1);
    }
}

int32_t XRand::NextIntUnbiased(int32_t a, int32_t b) {
    assert(a <= b);
    return a + NextUIntUnbiased(b - a + 1);
}

int64_t XRand::NextLongUnbiased(int64_t a, int64_t b) {
    assert(a <= b);
    return a + NextULongUnbiased(b - a + 1);
}

double XRand::NextDouble() {
    static const uint64_t kNumValues = 1LL << 53;
    return static_cast<double>(NextULong() & (kNumValues - 1)) / kNumValues;
}

double XRand::NextGaussian() {
    static const double kPi = acos(-1.0);
    const double value1 = NextDouble();
    const double value2 = NextDouble();
    return sqrt(-2.0 * log1p(-value1)) * cos(2.0 * kPi * value2);
}

template<typename RandomAccessIterator>
void XRand::Shuffle(RandomAccessIterator first, RandomAccessIterator last) {
    for (RandomAccessIterator iter = first; iter != last; ++iter) {
        std::iter_swap(first + NextUInt(iter - first + 1), iter);
    }
}

#endif  // #ifndef XRAND_H_

// 末尾に使用例
// 好きに改造して使ってください

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

class Graphics {
public:
    double screenW;
    double screenH;
    ostringstream data;

    double sr;
    double sg;
    double sb;
    double sa;
    double fr;
    double fg;
    double fb;
    double fa;

    Graphics() : screenW(1), screenH(1), sr(0), sg(0), sb(0), sa(1), fr(1), fg(1), fb(1), fa(1) {
    }

    tuple<double, double, double> color(double val) {
        if (val < 0.5) {
            double x = val * 2.0;
            return {30. * (1.0 - x) + 144. * x, 144. * (1.0 - x) + 255. * x, 255. * (1.0 - x) + 30. * x};
        } else {
            double x = val * 2.0 - 1.0;
            return {144. * (1.0 - x) + 255. * x, 255. * (1.0 - x) + 30. * x, 30. * (1.0 - x) + 70. * x};
        }
    }

    void screen(int width, int height) {
        screenW = width;
        screenH = height;
    }

    void clear() {
        data.str("");
        data.clear(stringstream::goodbit);
    }

    void stroke(double r, double g, double b) {
        stroke(r, g, b, 1);
    }

    void stroke(double r, double g, double b, double a) {
        sr = r;
        sg = g;
        sb = b;
        sa = a;
    }

    void noStroke() {
        stroke(0, 0, 0, 0);
    }

    void fill(double r, double g, double b) {
        fill(r, g, b, 1);
    }

    void fill(double r, double g, double b, double a) {
        fr = r;
        fg = g;
        fb = b;
        fa = a;
    }

    void noFill() {
        fill(0, 0, 0, 0);
    }

    void line(double x1, double y1, double x2, double y2) {
        data << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" " << stroke()
             << "/>\n";
    }

    void rect(double x, double y, double w, double h) {
        data << "<rect x=\"" << x << "\" y=\"" << y << "\" width=\"" << w << "\" height=\"" << h << "\" " << stroke()
             << " " + fill() << "/>\n";
    }

    void circle(double cx, double cy, double r) {
        data << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << r << "\" " << stroke() << " " + fill()
             << "/>\n";
    }

    void grid(vector <vector<int>> grid) {
        int h = grid.size();
        int w = grid[0].size();
        const int len = 10;
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                // 枠を黒に、塗りをいい感じに
                stroke(0.5, 0.5, 0.5);
                auto [r, g, b] = color((double) grid[i][j] / 100.0);
                fill(r / 256.0, g / 256.0, b / 256.0);
                rect(j * len, i * len + len / 2, len, len);
                // circle(j * len + len / 2, i * len + len, len / 2);
                // 文字を記述
                fill(0, 0, 0);
                text(to_string(grid[i][j]), j * len + len / 2, i * len + 2 * len / 3 + len / 2, 4);
            }
        }
    }

    void text(string str, double x, double y, double size = 16) {
        data << "<text text-anchor=\"middle\" x=\"" << x << "\" y=\"" << y << "\" font-size=\"" << size << "\" "
             << fill() << " >" << str << "</text>\n";
    }

    string dump(string id = "", string style = "", int widthPx = -1, int heightPx = -1) const {
        ostringstream res;
        res << "<svg ";
        if (id != "") res << "id=\"" + id + "\" ";
        if (style != "") res << "style=\"" + style + "\" ";
        if (widthPx != -1) res << "width=\"" << widthPx << "\" ";
        if (heightPx != -1) res << "height=\"" << heightPx << "\" ";
        res << "viewBox=\"-1 -1 " << screenW + 2 << " " << screenH + 2 << "\" xmlns=\"http://www.w3.org/2000/svg\">\n"
            << data.str() << "</svg>";
        return res.str();
    }

private:
    string stroke() const {
        return "stroke=\"" + rgb(sr, sg, sb) + "\" stroke-opacity=\"" + s(sa) + "\"";
    }

    string fill() const {
        return "fill=\"" + rgb(fr, fg, fb) + "\" fill-opacity=\"" + s(fa) + "\"";
    }

    string rgb(double r, double g, double b) const {
        return "rgb(" + s(lround(r * 255)) + "," + s(lround(g * 255)) + "," + s(lround(b * 255)) + ")";
    }

    string s(double x) const {
        return to_string(x);
    }
};

class Movie {
public:
    vector <string> svgs;

    Movie() {
    }

    void clear() {
        svgs.clear();
    }

    void addFrame(Graphics &g) {
        svgs.push_back(g.dump("f" + to_string(svgs.size()), "display:none;pointer-events:none;user-select:none;"));
    }

    string dumpHtml(int fps) {
        ostringstream res;
        res << "<html><body><div id=\"text\">loading...</div>" << endl;
        for (string &svg: svgs) {
            res << svg << endl;
        }

        res << "<script>\nlet numFrames = " << svgs.size() << ", fps = " << fps << ";";
        res << R"(
	let text = document.getElementById("text");
	let frames = [];
	for (let i = 0; i < numFrames; i++) {
		let f = document.getElementById("f" + i);
		frames.push(f);
		f.style.display = "none";
	}
	let currentFrame = 0;
	let playing = true;
	setInterval(() => {
		if (!playing) return;
		text.innerText = (currentFrame + 1) + " / " + numFrames;
		frames[(currentFrame - 1 + numFrames) % numFrames].style.display = "none";
		frames[currentFrame].style.display = null;
		currentFrame = (currentFrame + 1) % numFrames;
		if (currentFrame == 0) playing = false;
	}, 1000 / fps);
	window.onmousedown = e => { if (e.button == 0) playing = true; };
;)";
        res << "</script>" << endl;
        res << "</body></html>" << endl;
        return res.str();
    }

private:
};

int main() {
    Graphics g;
    Movie mov;
    int w = 500;
    int h = 500;

    XRand rnd(283);

    // スクリーンの大きさを設定
    g.screen(w, h);

    int n = 20;
    vector <vector<int>> grid(n, vector<int>(n, 0));

    for (int i = 0; i < 60 * 5; i++) {
        // 画面消去
        g.clear();
        // グリッドを描画
        for (int j = 0; j < 10; j++) {
            grid[rnd.NextInt(0, n - 1)][rnd.NextInt(0, n - 1)] += rnd.NextInt(1, 10);
        }
        g.grid(grid);
        mov.addFrame(g);
    }
    // 60FPS で書き出し
    string html = mov.dumpHtml(60);

    // 動画を保存
    ofstream fout;
    fout.open("movie.html", ios::out);
    if (!fout) {
        return 1;
    }
    fout << html << endl;
    fout.close();

    return 0;
}
