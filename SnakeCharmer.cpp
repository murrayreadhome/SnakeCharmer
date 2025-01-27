// C++11
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <random>
#include <cmath>
#include <numeric>
#include <deque>
#include <queue>
#include <list>
#include <complex>
#include <iomanip>
#include <iterator>
#include <cstring>
#include <memory>

using namespace std;

#ifdef TESTING_AT_HOME
#ifdef _DEBUG
#define TIME_LIMIT 100000
#define MIN_STEPS 0
#else
#define TIME_LIMIT 8000
#define MIN_STEPS 0
#endif
#else
#define TIME_LIMIT 8000
#define MIN_STEPS 0
#endif

#ifdef TESTING_AT_HOME
#include <time.h>
int milliseconds(bool reset = false)
{
    static clock_t start = clock();
    if (reset)
        start = clock();
    clock_t now = clock();
    clock_t elapsed = now - start;
    elapsed = (elapsed * 1000) / (CLOCKS_PER_SEC);
    return int(elapsed);
}
#else
#include <sys/time.h>
int getTime()
{
    timeval t;
    gettimeofday(&t, NULL);
    int tm = t.tv_sec;
    tm *= 1000;
    tm += (t.tv_usec / 1000);
    return tm;
}
int milliseconds(bool /*reset*/ = false)
{
    static int start = getTime();
    int now = getTime();
    if (now < start)
        now += 60000;	// to account for minute rollover
    int elapsed = now - start;
    return elapsed;
}
#endif

#define ALL(cont) (cont).begin(), (cont).end()

template<bool Condition> struct STATIC_ASSERTION_FAILURE { enum { ASize = -1 }; };
template <> struct STATIC_ASSERTION_FAILURE<true> { enum { ASize = 1 }; };
#define STATIC_ASSERT(cond, name) extern int name[STATIC_ASSERTION_FAILURE<cond>::ASize];

bool keepGoing(int aStart, int aLimit, int aSteps, int aMaxSteps = 2000000000, int aMinSteps = MIN_STEPS)
{
    // max takes precedence over min
    if (aSteps >= aMaxSteps)
        return false;
    if (aSteps < aMinSteps)
        return true;
    int now = milliseconds();
    if (now > aLimit)
        return false;
    int perStep = (now - aStart) / max(aSteps, 1);
    //	cerr << "perstep = " << perStep << " est = " << (aStart + perStep * (aSteps + 1)) << endl;
    return (aStart + perStep * (aSteps + 1)) < aLimit;
}

typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char uchar;
typedef long long LL;

struct Error : public std::exception
{
    string iErr;
    Error(const string& e)
        : iErr(e)
    {}
    ~Error() throw()
    {}
    const char* what() const throw()
    {
        return iErr.c_str();
    }
};

void report(const string& msg)
{
#ifdef TESTING_AT_HOME
    throw Error(msg);
#else
    cerr << msg << endl;
#endif
}

template <typename T>
bool incr(const T& a, const T& b, const T& c)
{
    return a <= b && b <= c;
}
template <typename T>
bool incrStrict(const T& a, const T& b, const T& c)
{
    return a < b && b < c;
}
template <typename T>
bool between(const T& a, const T& b, const T& c)
{
    return min(a, c) <= b && b <= max(a, c);
}
template <typename T>
bool betweenStrict(const T& a, const T& b, const T& c)
{
    return min(a, c) < b && b < max(a, c);
}

const int KUp = -1;
enum TDir { EUp = 0, ERight, EDown, ELeft, ENumDirs, ENorth = EUp, EEast = ERight, ESouth = EDown, EWest = ELeft, ENoDir = 8 };
bool Horizontal(TDir d) { return d == ERight || d == ELeft; }

template <typename T>
struct Point2D
{
    T x;
    T y;
    Point2D()
        : x(), y()
    {}
    template <typename U1, typename U2>
    Point2D(const U1& aX, const U2& aY)
        : x(aX), y(aY)
    {}
    template <typename U>
    Point2D(const Point2D<U> aPoint2D)
        : x(T(aPoint2D.x)), y(T(aPoint2D.y))
    {}
    template <typename U>
    const Point2D& operator=(const Point2D<U> aPoint2D)
    {
        x = T(aPoint2D.x);
        y = T(aPoint2D.y);
        return *this;
    }
    template <typename U>
    Point2D operator+(const Point2D<U>& aPoint2D) const
    {
        return Point2D(x + aPoint2D.x, y + aPoint2D.y);
    }
    template <typename U>
    const Point2D& operator+=(const Point2D<U>& aPoint2D)
    {
        x += aPoint2D.x;
        y += aPoint2D.y;
        return *this;
    }
    Point2D operator-(const Point2D& aPoint2D) const
    {
        return Point2D(x - aPoint2D.x, y - aPoint2D.y);
    }
    const Point2D& operator-=(const Point2D& aPoint2D)
    {
        x -= aPoint2D.x;
        y -= aPoint2D.y;
        return *this;
    }
    bool operator==(const Point2D& aPoint2D) const
    {
        return x == aPoint2D.x && y == aPoint2D.y;
    }
    bool operator!=(const Point2D& aPoint2D) const
    {
        return x != aPoint2D.x || y != aPoint2D.y;
    }
    Point2D operator*(T aFactor) const
    {
        return Point2D(x * aFactor, y * aFactor);
    }
    const Point2D& operator*=(T aFactor)
    {
        x = x * aFactor;
        y = y * aFactor;
        return *this;
    }
    Point2D operator/(T aFactor) const
    {
        return Point2D(x / aFactor, y / aFactor);
    }
    const Point2D& operator/=(T aFactor)
    {
        x = x / aFactor;
        y = y / aFactor;
        return *this;
    }
    Point2D operator-() const
    {
        return Point2D(-x, -y);
    }
    bool operator<(const Point2D& aOther) const
    {
        return (y < aOther.y) || ((y == aOther.y) && (x < aOther.x));
    }
    T cross(const Point2D& aP) const
    {
        return x * aP.y - y * aP.x;
    }
    T operator*(const Point2D& aP) const
    {
        return cross(aP);
    }
    T dot(const Point2D& aP) const
    {
        return x * aP.x + y * aP.y;
    }
    T operator()(const Point2D& aP) const
    {
        return dot(aP);
    }
    T boxDist(const Point2D& aPoint2D) const
    {
        return abs(x - aPoint2D.x) + abs(y - aPoint2D.y);
    }
    T squaredMagnitude() const
    {
        return x * x + y * y;
    }
    double magnitude() const
    {
        return sqrt(double(squaredMagnitude()));
    }
    double dist(const Point2D& aPoint2D) const
    {
        T dx = x - aPoint2D.x;
        T dy = y - aPoint2D.y;
        return sqrt(double(dx * dx + dy * dy));
    }
    double angle()
    {
        return atan2(double(y), double(x));
    }
    Point2D<double> unit() const
    {
        double mag = magnitude();
        if (mag != 0.0)
            return Point2D<double>(x / mag, y / mag);
        else
            return Point2D<double>(0.0, 0.0);
    }
    Point2D rotate(double aAngle) const
    {
        double c = cos(aAngle);
        double s = sin(aAngle);
        return Point2D(T(x * c - y * s), T(x * s + y * c));
    }
    Point2D rotateCW() const
    {
        return Point2D(-y, x);
    }
    Point2D rotateCCW() const
    {
        return Point2D(y, -x);
    }
    Point2D swapped() const
    {
        return Point2D(y, x);
    }
    T operator[](TDir dir) const
    {
        switch (dir)
        {
        case EUp:
        case EDown:
            return y;
        case ELeft:
        case ERight:
            return x;
        default:
            return T();
        }
    }
    T& operator[](TDir dir)
    {
        switch (dir)
        {
        case EUp:
        case EDown:
            return y;
        case ELeft:
        case ERight:
        default:
            return x;
        }
    }
};

template <typename T>
Point2D<T> minPos(const Point2D<T>& a, const Point2D<T>& b)
{
    return Point2D<T>(min(a.x, b.x), min(a.y, b.y));
}

template <typename T>
Point2D<T> maxPos(const Point2D<T>& a, const Point2D<T>& b)
{
    return Point2D<T>(max(a.x, b.x), max(a.y, b.y));
}

template <typename T>
ostream& operator<<(ostream& aStream, const Point2D<T>& aPoint2D)
{
    aStream << aPoint2D.x << "," << aPoint2D.y;
    return aStream;
}

typedef Point2D<int> Pos;
typedef Point2D<LL> PosLL;
typedef Point2D<double> Point;
typedef Point2D<float> PointF;
template <typename T>
Pos toPos(const Point2D<T>& p)
{
    return Pos(int(p.x + 0.5), int(p.y + 0.5));
}

const double triEps = 1e-12;
template <typename T>
bool pointInTriangle(const Point2D<T>& p, const Point2D<T>& a, const Point2D<T>& b, const Point2D<T>& c)
{
    double x = (a - b) * (a - p), y = (b - c) * (b - p), z = (c - a) * (c - p);
    return ((x <= triEps && y <= triEps && z <= triEps) || (x >= -triEps && y >= -triEps && z >= -triEps));
}

const Pos KDir[9] = { Pos(0,KUp), Pos(1,0), Pos(0,-KUp), Pos(-1,0), Pos(1,-KUp), Pos(1,KUp), Pos(-1,-KUp), Pos(-1,KUp), Pos(0,0) };
const string KDirName = "URDL";
const string KCompassName = "NESW";
const TDir KOppositeDir[4] = { EDown, ELeft, EUp, ERight };

TDir dirTo(const Pos& aFrom, const Pos& aTo)
{
    int dx = aTo.x - aFrom.x;
    int dy = aTo.y - aFrom.y;
    //	if (dx && dy)
    //		throw Error("No Direction!");
    if (abs(dy) > abs(dx))
    {
        bool b = dy < 0;
#ifdef TESTING_AT_HOME
#pragma warning(disable:4127)
#endif
        if (KUp < 0)
            b = !b;
#ifdef TESTING_AT_HOME
#pragma warning(default:4127)
#endif
        return b ? EDown : EUp;
    }
    else
        return dx > 0 ? ERight : (dx < 0 ? ELeft : ENoDir);
}

struct TDirs { TDir x, y; };
TDirs dirsTo(const Pos& aFrom, const Pos& aTo)
{
    TDirs dirs = { ENoDir, ENoDir };
    if (aTo.x < aFrom.x)
        dirs.x = ELeft;
    if (aTo.x > aFrom.x)
        dirs.x = ERight;
    if (aTo.y < aFrom.y)
        dirs.y = EUp;
    if (aTo.y > aFrom.y)
        dirs.y = EDown;
#if	KUp == 1
    dirs.second = TDir(EDown - dirs.second);
#endif
    return dirs;
}

template <typename T>
struct Grid
{
    typedef vector<T> Coll;
    typedef typename Coll::iterator iterator;
    typedef typename Coll::const_iterator const_iterator;
    typedef typename Coll::reference Ref;
    typedef typename Coll::const_reference ConstRef;

    Grid()
        : width(0), height(0), bWidth(0), bHeight(0), border(0)
    {}

    Grid(int aWidth, int aHeight, const T& aVal)
    {
        init(aWidth, aHeight, aVal, 0);
    }

    Grid(int aWidth, int aHeight, const T& aVal, int aBorder)
    {
        init(aWidth, aHeight, aVal, aBorder);
    }

    Grid(const Grid& aGrid)
    {
        *this = aGrid;
    }

    const Grid& operator=(const Grid& aGrid)
    {
        if (&aGrid == this)
            return *this;
        width = aGrid.width;
        height = aGrid.height;
        bWidth = aGrid.bWidth;
        bHeight = aGrid.bHeight;
        border = aGrid.border;
        iGrid = aGrid.iGrid;
        iBegin = iGrid.begin();
        iBegin = iter(Pos(border, border));
        return *this;
    }

    void init(int aWidth, int aHeight, const T& aVal, int aBorder = 0)
    {
        bWidth = aWidth + aBorder * 2;
        bHeight = aHeight + aBorder * 2;
        width = aWidth;
        height = aHeight;
        border = aBorder;
        iGrid.clear();
        iGrid.resize(bWidth * bHeight, aVal);
        iBegin = iGrid.begin();
        iBegin = iter(Pos(aBorder, aBorder));
    }

    inline int index(const Pos& aPos) const
    {
#ifdef _DEBUG
        size_t base = border * (bWidth + 1);
        size_t idx = aPos.x + aPos.y * bWidth + base;
        if (idx < 0 || idx > size(iGrid))	// allow for "end" positions
            report("Grid bad pos");
        return int(idx - base);
#else
        return aPos.x + aPos.y * bWidth;
#endif
    }

    Ref operator[](const Pos& aPos)
    {
        return iBegin[index(aPos)];
    }

    T operator[](const Pos& aPos) const
    {
        return iBegin[index(aPos)];
    }

    ConstRef get(const Pos& aPos) const
    {
        return iBegin[index(aPos)];
    }

    Ref get(const Pos& aPos)
    {
        return iBegin[index(aPos)];
    }

    void set(const Pos& aPos, const T& v)
    {
        iBegin[index(aPos)] = v;
    }

    iterator begin()
    {
        return iGrid.begin();
    }

    iterator end()
    {
        return iGrid.end();
    }

    const_iterator begin() const
    {
        return iGrid.begin();
    }

    const_iterator end() const
    {
        return iGrid.end();
    }

    iterator iter(const Pos& aPos)
    {
        return iBegin + index(aPos);
    }

    const_iterator iter(const Pos& aPos) const
    {
        return iBegin + index(aPos);
    }

    bool isValid(const Pos& aPos) const
    {
        return 0 <= aPos.x && aPos.x < width && 0 <= aPos.y && aPos.y < height;
    }

    bool isBorderValid(const Pos& aPos) const
    {
        return -border <= aPos.x && aPos.x < width + border && -border <= aPos.y && aPos.y < height + border;
    }

    void limit(Pos& aPos) const
    {
        aPos.x = min(max(0, aPos.x), width - 1);
        aPos.y = min(max(0, aPos.y), height - 1);
    }

    template<typename T2>
    void swap(Grid<T2>& other)
    {
        iGrid.swap(other.iGrid);
        swap(width, other.width);
        swap(height, other.height);
        swap(bWidth, other.bWidth);
        swap(bHeight, other.bHeight);
        swap(border, other.border);
        swap(iBegin, other.iBegin);
    }

    int width;
    int height;
    int bWidth;
    int bHeight;
    int border;
    Coll iGrid;
    typename Coll::iterator iBegin;
};

#define FOR_GRID(var, grid)    for (Pos var(0,0); var.y<grid.height; var.y++) for (var.x=0; var.x<grid.width; var.x++)
#define FOR_SIZE(var, w, h)    for (Pos var(0,0); var.y<(h); var.y++) for (var.x=0; var.x<(w); var.x++)
#define FOR_RECT(var, l, t, r, b)    for (Pos var((l),(t)); var.y<(b); var.y++) for (var.x=(l); var.x<(r); var.x++)

template<typename T>
ostream& operator<<(ostream& aStream, const Grid<T>& aGrid)
{
    for (Pos var(0, 0); var.y < aGrid.height; var.y++)
    {
        for (var.x = 0; var.x < aGrid.width; var.x++)
        {
            aStream << aGrid[var] << ",";
        }
        aStream << endl;
    }
    return aStream;
}

template<typename T>
int Size(const T& cont)
{
    return int(cont.size());
}

template<typename T>
ostream& operator<<(ostream& aStream, const vector<T>& aVect)
{
    int s = Size(aVect);
    aStream << "{";
    for (int i = 0; i < s; i++)
    {
        if (i > 0)
            aStream << ",";
        aStream << aVect[i];
    }
    aStream << "}" << endl;
    return aStream;
}

const double PI = acos(-1.0);
const double PI2 = PI * 2;
const double SQRT2 = sqrt(2.0);

inline double angleNorm(double a)
{
    if (a < -PI)
    {
        int times = int((a - PI) / PI2);
        a -= times * PI2;
    }
    else if (a > PI)
    {
        int times = int((a + PI) / PI2);
        a -= times * PI2;
    }
    return a;
}
inline double angleAdd(double a, double b)
{
    return angleNorm(a + b);
}
inline double angleSub(double a, double b)
{
    return angleNorm(a - b);
}

struct OrderPointByAngleToRef
{
    Point iRef;
    OrderPointByAngleToRef(const Point& aRef) : iRef(aRef)
    {}
    bool operator()(const Point& aL, const Point& aR) const
    {
        return (aL - iRef).angle() > (aR - iRef).angle();
    }
};

struct OrderPointByDistToRef
{
    Point iRef;
    OrderPointByDistToRef(const Point& aRef) : iRef(aRef)
    {}
    bool operator()(const Point& aL, const Point& aR) const
    {
        return (aL - iRef).squaredMagnitude() < (aR - iRef).squaredMagnitude();
    }
};

uint hash(const string& str)
{
    uint hash = 5381;
    for (int i = 0; i < Size(str); i++)
        hash = hash * 33 ^ str[i];
    return hash;
}

int bitCount[256] =
{
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

inline int countBits(int bits)
{
    unsigned char* bytes = (unsigned char*)&bits;
    return bitCount[bytes[0]] + bitCount[bytes[1]] + bitCount[bytes[2]] + bitCount[bytes[3]];
}

template<typename TDATA, typename TSCORE>
class TrackMin
{
public:
    TDATA data;
    TSCORE val;
    TrackMin(const TDATA& startData, const TSCORE& startVal)
        : data(startData), val(startVal)
    {}

    bool update(const TDATA& newData, const TSCORE& newVal)
    {
        if (newVal < val)
        {
            data = newData;
            val = newVal;
            return true;
        }
        return false;
    }
};

template<typename TDATA, typename TSCORE>
class TrackMax
{
public:
    TDATA data;
    TSCORE val;
    TrackMax(const TDATA& startData, const TSCORE& startVal)
        : data(startData), val(startVal)
    {}

    bool update(const TDATA& newData, const TSCORE& newVal)
    {
        if (newVal > val)
        {
            data = newData;
            val = newVal;
            return true;
        }
        return false;
    }
};

template<typename TDATA, typename TSCORE>
class TrackMaxEq
{
public:
    TDATA data;
    TSCORE val;
    TrackMaxEq(const TDATA& startData, const TSCORE& startVal)
        : data(startData), val(startVal)
    {}

    bool update(const TDATA& newData, const TSCORE& newVal)
    {
        if (newVal >= val)
        {
            data = newData;
            val = newVal;
            return true;
        }
        return false;
    }
};

template <typename T, typename U>
U& other(const T& from, U& a, U& b)
{
    if (from == a)
        return b;
    else if (from == b)
        return a;
    else
        report("neither");
    return a;
}

#ifdef TESTING_AT_HOME
inline bool isnan(double x)
{
    return x != x;
}
#endif

////////////////////////////////////////////////////////
////////////////  Solution starts here  ////////////////
////////////////////////////////////////////////////////

const static Pos char_dirs[96] = 
{
    {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
    {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
    {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
    {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},
    {},{},{},{},KDir[EDown],{},{},{},{},{},{},{},KDir[ELeft],{},{},{},
    {},{},KDir[ERight],{},{},KDir[EUp],{},{},{},{},{},{},{},{},{},{},
};


struct Problem
{
    int N;
    int V;
    string snake;

    void generate(int seed)
    {
        mt19937 re(seed);
        N=uniform_int_distribution<int>(3,24)(re)*2 + 1;
        V=uniform_int_distribution<int>(2,8)(re);
        snake.resize(N*N);
        for (size_t i=0; i<size_t(N*N); i++)
            snake[i] = '0' + uniform_int_distribution<int>(2,V+1)(re);
    }

    mutable Grid<int> g;
    int score_chars(const vector<char>& path, bool print_grid = false) const
    {
        Pos p(N / 2, N / 2);
        g.init(N, N, 0, 1);
        size_t i = path.size();
        g[p] = snake[i--]-'0';
        for (char c : path)
        {
            p += char_dirs[c];
            if (!g.isValid(p))
                report("Score out of grid");
            if (g[p] != 0)
                report("Score crash");
            g[p] = snake[i--] - '0';
        }

        int score = 0;
        FOR_GRID(p, g)
        {
            int cell = g[p];
            int prod = cell;
            for (int d = 0; d < 4; d++)
                if (g[p + KDir[d]] == cell)
                    prod *= cell;
            score += prod;
        }

        if (print_grid)
            cout << g;

        return score;
    }
};

istream& operator>>(istream& in, Problem& p)
{
    in >> p.N;
    in >> p.V;
    in >> p.snake;
    return in;
}

ostream& operator<<(ostream& out, const Problem& p)
{
    out << p.N << endl;
    out << p.V << endl;
    out << p.snake << endl;
    return out;
}

class PlacementFail : public std::exception 
{
public:
    size_t id;
    PlacementFail(size_t id=0) : id(id) {}
};

class SnakeCharmer 
{
    size_t N;
    int maxV;
    vector<int> snake;
    size_t n;
    Problem problem;

public:
    vector<char> findSolution(int aN, int aV, string aSnake)
    {
        milliseconds();

        init(aN, aV, aSnake);
        simple_solution();
        find_best_snake();
        
        return best_path;
    }

    void init(int aN, int aV, string aSnake)
    {
        problem = Problem{ aN, aV, aSnake };
        N = aN;
        maxV = aV + 1;
        snake.reserve(aSnake.size());
        for (char c : aSnake)
            snake.push_back(c - '0');
        n = snake.size();
    }

    void simple_solution()
    {
        char names[] = { 'L','D','R','U' };

        int n = int(N * N - 1);
        vector<char> moves(n);

        for (int i = 0, dir = 0, L = 1; i < n; )
        {
            for (int k = 0; k < L && i < n; k++, i++) moves[i] = names[dir];
            dir = (dir + 1) % 4;
            for (int k = 0; k < L && i < n; k++, i++) moves[i] = names[dir];
            dir = (dir + 1) % 4;
            L++;
        }

        best_path = moves;
        best_score = problem.score_chars(moves);
    }

    int bestV;
    size_t bestN;
    int bestC;
    void find_best_snake()
    {
        for (int v = maxV; v > (maxV+1)/2; v--)
        {
            size_t last_start = n;
            find_runs(v);
            int reductions = 0;
            for (size_t n = N - 1; n > 0; n--)
            {
                auto ms = milliseconds();
                //cout << ms << endl;
                if (ms > TIME_LIMIT)
                    break;
                vector<pair<size_t, size_t>> meeting = best_meeting(n, v, last_start);
                if (meeting.size() < n)
                    n = meeting.size();
                if (n == 0)
                    break;
                reductions++;
                if (reductions > 3)
                    last_start = non_v_before(meeting.back().first, v);
                for (int c = 7; c >= 0; c--)
                {
                    // cout << v << " " << n << " " << c << endl;
                    try
                    {
                        initial_placement(meeting, c, v);
                        //cout << "a" << endl;
                        connect_snake();
                        size_t score = best_score;
                        score_and_record();
                        if (score != best_score)
                        {
                            bestV = v;
                            bestN = n;
                            bestC = c;
                        }
                    }
                    catch (PlacementFail)
                    {
                    //    cout << "fail" << endl;
                    }
                }
            }
        }
    }

    struct Run
    {
        size_t start;
        size_t length;
    };
    vector<Run> runs;

    void find_runs(int v)
    {
        runs.clear();

        size_t run_start = n;
        auto record_run = [&](size_t i) 
        { 
            runs.push_back({ run_start, i - run_start }); 
        };

        for (size_t i = 0; i < n; i++)
        {
            if (snake[i] == v)
            {
                if (run_start == n)
                    run_start = i;
            }
            else if (run_start != n)
            {
                record_run(i);
                run_start = n;
            }
        }

        if (run_start != n)
            record_run(n);
    }

    size_t score_runs(const vector<Run>& runs, size_t v)
    {
        size_t total = 0;
        for (size_t i=0; i<runs.size(); i++)
        {
            Run a{0,0};
            if (i>0)
                a=runs[i-1];
            Run b=runs[i];
            Run c{0,0};
            if (i < runs.size()-1)
                c=runs[i+1];

            size_t parity = i%2;

            size_t a_offset = (a.start%2 != parity) ? 0 : 1;
            size_t b_offset = (b.start%2 == parity) ? 0 : 1;
            size_t c_offset = (c.start%2 != parity) ? 0 : 1;

            for (size_t j=0; j<b.length; j++)
            {
                size_t cell = v;
                if (j>0) cell *= v;
                if (j<b.length-1) cell *= v;
                size_t a_rel = j + a.length/2 + a_offset - b.length/2 - b_offset;
                if (0 <= a_rel && a_rel < a.length) cell *= v;
                size_t c_rel = j + c.length/2 + c_offset - b.length/2 - b_offset;
                if (0 <= c_rel && a_rel < c.length) cell *= v;
                total += cell;
            }
        }
        return total;
    }

    vector<pair<size_t, size_t>> best_meeting(size_t n, int v, size_t max_start)
    {
        vector<pair<size_t, size_t>> best;
        size_t best_score = 0;
        for (size_t parity = 0; parity < 2; parity++)
        {
            vector<Run> valid_runs;
            for (const Run& run : runs)
            {
                if (run.start >= max_start)
                    break;
                if (run.length == 2)
                {
                    if ((run.start % 2) == parity)
                        valid_runs.push_back(run);
                }
                else if (run.length > 2)
                {
                    size_t start = run.start + run.length / 2 - 1;
                    if ((start % 2) == parity)
                        valid_runs.push_back({ start, run.length });
                    else
                        valid_runs.push_back({ start+1, run.length });  // todo extra test for optimal 4+ length offset of +/-1
                }
            }

            while (valid_runs.size() > n)
            {                
                // todo could attempt to bring large runs together
                // removing from the end, to free up more space, but does that really help? or just cost a little bit?
                auto smallest = min_element(valid_runs.rbegin(), valid_runs.rend(), [](const Run& a, const Run& b)
                {
                    return a.length < b.length;
                });
                if (smallest->length > 2)
                    valid_runs.resize(n);
                else
                    valid_runs.erase((++smallest).base());
            }

            // todo improve scoring, size probably matters, but isn't everything
            size_t runs_score = 0;
            for (const Run& run : valid_runs)
                runs_score += run.length;

//            size_t runs_score = score_runs(valid_runs, v);
            if (runs_score > best_score)
            {
                best_score = runs_score;
                best.clear();
                for (const Run& run : valid_runs)
                    best.push_back({ run.start, run.start + 1 });
            }
        }

        return best;
    }

    Grid<size_t> placement;

    struct Growth
    {
        size_t where;
        vector<Pos> add;
    };

    struct Loop
    {
        size_t id;
        size_t first;
        size_t last;
        vector<Pos> points;
        bool start_free;
        bool single_growth;
        //Pos target;
        size_t size() const { return 1 + last - first; }
        size_t remaining() const { return size() - points.size(); }
        bool complete() const { return size() == points.size(); }
        Grid<double> pressure;
        Grid<int> visited;
        vector<Growth> growth;
    };
    vector<Loop> loops;

    void initial_placement(const vector<pair<size_t, size_t>>& meeting, int cap_options, int v)
    {
        if (meeting.empty())
            throw PlacementFail();

        placement.init(N, N, n, 1);
        loops.clear();

        int left_size_adjust = 0;
        if (cap_options == 7)
        {
            left_size_adjust = -1;
            cap_options = 6;
        }
        bool end_caps = cap_options & 1;
        bool mid_caps = cap_options & 2;
        bool odd_caps = cap_options & 4;

        size_t parity = meeting.front().first % 2;
        size_t m = meeting.size();
        size_t m2 = m / 2;

        Pos p(N / 2, N / 2);
        size_t loop_end = n;
        Pos loop_end_pos = p;
        size_t meeting_idx = m - 1;
        bool ignore_next_loop = false;

        auto place_connection = [&](TDir d1, TDir d2)
        {
            const auto& connection = meeting[meeting_idx--];
            if (!placement.isValid(p))
                throw PlacementFail();
            placement[p] = connection.second;
            if (loop_end != n && !ignore_next_loop)
            {
                Loop loop = { loops.size(), connection.second, loop_end, {p, loop_end_pos}, false };
                loops.push_back(loop);
            }
            ignore_next_loop = false;
            p += KDir[d1];
            if (!placement.isValid(p))
                throw PlacementFail();
            placement[p] = connection.first;
            loop_end = connection.first;
            loop_end_pos = p;
            p += KDir[d2];
        };

        auto place_cap = [&](size_t cap, TDir dir)
        {
            if (cap != n)
            {
                placement[p] = cap;
                if (loop_end != n)
                {
                    Loop loop = { loops.size(), cap, loop_end, {p, loop_end_pos}, false };
                    loops.push_back(loop);
                }
                loop_end = cap;
                loop_end_pos = p;
            }
            p += KDir[dir];
        };

        auto add_mid_caps = [&](TDir d1, TDir d2)
        {
            if (!placement.isValid(p))
                throw PlacementFail();
            size_t lower_end = meeting[meeting_idx].second;
            size_t gap_start = non_v_after(lower_end + 1, v);
            size_t gap_end = non_v_before(loop_end - 1, v);
            size_t low_cap = find_between(gap_start, gap_end, v, parity);
            size_t high_cap = find_between(low_cap == n ? gap_start : low_cap + 1, gap_end, v, 1 - parity);
            place_cap(high_cap, d1);
            place_cap(low_cap, d2);
        };

        auto odd_cap_connection = [&](TDir d1, TDir d2)
        {
            Pos cap1pos = p - KDir[d2];
            Pos join1start = cap1pos - KDir[d2];
            if (!placement.isValid(join1start))
                throw PlacementFail();
            const auto& connection = meeting[meeting_idx];
            size_t cap1 = find_between(non_v_after(connection.second + 1, v), non_v_before(loop_end - 1, v), v, parity);
            size_t inner_len = cap1 - connection.second;
            size_t outer_len = loop_end - cap1;
            if (inner_len+5 > outer_len)
            {
                // cap loops won't fit
                place_connection(d1, d2);
            }
            else
            {
                size_t fold_len = inner_len / 2;
                Loop outerloop = { loops.size(), cap1, loop_end, {cap1pos}, false };
                Pos lp = join1start;
                outerloop.points.push_back(lp);
                placement[lp] = outerloop.id;
                for (size_t i = 0; i < fold_len + 1; i++)
                {
                    lp -= KDir[d1];
                    if (!placement.isValid(lp))
                        throw PlacementFail();
                    outerloop.points.push_back(lp);
                    placement[lp] = outerloop.id;
                }
                for (size_t i = 0; i < 3; i++)
                {
                    lp += KDir[d2];
                    outerloop.points.push_back(lp);
                    placement[lp] = outerloop.id;
                }
                for (size_t i = 0; i < fold_len + 1; i++)
                {
                    lp += KDir[d1];
                    outerloop.points.push_back(lp);
                    placement[lp] = outerloop.id;
                }
                loops.push_back(outerloop);

                Loop innerloop = { loops.size(), connection.second, cap1, {p}, false };
                lp = p;
                for (size_t i = 0; i < fold_len; i++)
                {
                    lp -= KDir[d1];
                    innerloop.points.push_back(lp);
                    placement[lp] = innerloop.id;
                }
                lp -= KDir[d2];
                innerloop.points.push_back(lp);
                placement[lp] = innerloop.id;
                for (size_t i = 0; i < fold_len; i++)
                {
                    lp += KDir[d1];
                    innerloop.points.push_back(lp);
                    placement[lp] = innerloop.id;
                }
                loops.push_back(innerloop);

                ignore_next_loop = true;
                place_connection(d1, d2);
            }

            if (meeting_idx < meeting.size())
            {
                const auto& next_connection = meeting[meeting_idx];
                Pos cap2pos = cap1pos + KDir[d1];
                Pos join2start = p;
                size_t cap2 = find_between_down(non_v_after(next_connection.second + 1, v), non_v_before(connection.first - 1, v), v, 1-parity);
                size_t inner_len = connection.first - cap2;
                size_t outer_len = cap2 - next_connection.second;
                if (inner_len + 5 <= outer_len)
                {
                    size_t fold_len = inner_len / 2;
                    Loop outerloop = { loops.size(), next_connection.second, cap2, {join2start}, false };
                    Pos lp = join2start;
                    for (size_t i = 0; i < fold_len + 1; i++)
                    {
                        lp += KDir[d1];
                        if (!placement.isValid(lp))
                            throw PlacementFail();
                        outerloop.points.push_back(lp);
                        placement[lp] = outerloop.id;
                    }
                    for (size_t i = 0; i < 3; i++)
                    {
                        lp -= KDir[d2];
                        outerloop.points.push_back(lp);
                        placement[lp] = outerloop.id;
                    }
                    for (size_t i = 0; i < fold_len + 1; i++)
                    {
                        lp -= KDir[d1];
                        outerloop.points.push_back(lp);
                        placement[lp] = outerloop.id;
                    }
                    loops.push_back(outerloop);

                    Loop innerloop = { loops.size(), cap2, connection.first, {cap2pos}, false };
                    lp = cap2pos;
                    placement[lp] = innerloop.id;
                    for (size_t i = 0; i < fold_len; i++)
                    {
                        lp += KDir[d1];
                        innerloop.points.push_back(lp);
                        placement[lp] = innerloop.id;
                    }
                    lp += KDir[d2];
                    innerloop.points.push_back(lp);
                    placement[lp] = innerloop.id;
                    for (size_t i = 0; i < fold_len; i++)
                    {
                        lp -= KDir[d1];
                        innerloop.points.push_back(lp);
                        placement[lp] = innerloop.id;
                    }
                    loops.push_back(innerloop);

                    ignore_next_loop = true;
                }
            }
        };

        if (end_caps)
        {
            size_t e_cap = find_between(meeting.back().second + 1, n - 1, v, parity);
            place_cap(e_cap, ELeft);

            for (size_t i = 0; i < m2; i++)
                place_connection(ELeft, ELeft);

            // odd middle goes vertical
            if (m % 2 == 1)
            {
                if (odd_caps)
                    odd_cap_connection(EUp, ERight);
                else
                    place_connection(EUp, ERight);
            }
            else if (mid_caps)
                add_mid_caps(EUp, ERight);
            else
                p = p + KDir[EUp] + KDir[ERight];

            for (size_t i = 0; i < m2; i++)
                place_connection(ERight, ERight);

            size_t s_cap = find_between(0, meeting.front().first - 1, v, 1-parity);
            place_cap(s_cap, ERight);
        }
        else
        {
            if (m < 3)
                throw PlacementFail();

            size_t left_size = 1;
            while (meeting.back().second / 2 < meeting[m - 1 - left_size].second)
                left_size++;
            left_size += left_size_adjust;
            size_t right_size = m - left_size;

            size_t m4 = left_size / 2;
            for (size_t i = 0; i < m4; i++)
                place_connection(ELeft, ELeft);

            if (mid_caps && left_size % 2 == 0)
                add_mid_caps(EUp, ERight);
            else if (left_size % 2 == 1)
            {
                if (odd_caps)
                    odd_cap_connection(EUp, ERight);
                else
                    place_connection(EUp, ERight);
            }
            else
                p = p + KDir[EUp] + KDir[ERight];

            for (size_t i = 0; i < m4; i++)
                place_connection(ERight, ERight);

            m4 = right_size / 2;
            for (size_t i = 0; i < m4; i++)
                place_connection(ERight, ERight);

            if (mid_caps && right_size % 2 == 0)
                add_mid_caps(EDown, ELeft);
            else if (right_size % 2 == 1)
            {
                if (odd_caps)
                    odd_cap_connection(EDown, ELeft);
                else
                    place_connection(EDown, ELeft);
            }
            else
                p = p + KDir[EDown] + KDir[ELeft];

            for (size_t i = 0; i < m4; i++)
                place_connection(ELeft, ELeft);
        }

        // add loop to start
        Loop loop = { loops.size(), 0, loop_end, {loop_end_pos}, true };
        loops.push_back(loop);
    }

    size_t non_v_after(size_t i, size_t v)
    {
        while (i<n && snake[i] == v)
            i++;
        return i;
    }

    size_t non_v_before(size_t i, size_t v)
    {
        while (i<n && snake[i] == v)
            i--;
        if (i > n)
            i = n;
        return i;
    }

    size_t find_between(size_t from, size_t to, int v, size_t parity)
    {
        for (size_t i = from; i <= to; i++)
        {
            if ((i % 2) == parity && snake[i] == v)
                return i;
        }
        return n;
    }

    size_t find_between_down(size_t from, size_t to, int v, size_t parity)
    {
        for (size_t i = to; n > i && i >= from; i--)
        {
            if ((i % 2) == parity && snake[i] == v)
                return i;
        }
        return n;
    }

    void connect_snake()
    {
        connect_cap_misses();
        do
        {
            analyse_loops();
            if (grow_single_options())
                continue;
            grow_highest_pressure();
        } while (!all_loops_done());
        place_loop_numbers();
    }

    void connect_cap_misses()
    {
        for (Loop& loop : loops)
        {
            if (loop.points.size() == 2 &&
                loop.points[0].boxDist(loop.points[1]) == 2)
            {
                Pos p(loop.points[0].x, loop.points[1].y);
                if (placement[p] != n)
                    p = Pos(loop.points[1].x, loop.points[0].y);
                if (placement[p] != n)
                    report("should not have a 2 gap");
                placement[p]=loop.id;
                loop.points.insert(loop.points.begin()+1, p);
            }
        }
    }

    void add_loop_growth_points(Loop& loop)
    {
        // find growth points
        for (size_t i = 1; loop.remaining() >= 2 && i < loop.points.size(); i++)
        {
            const Pos& a = loop.points[i - 1];
            const Pos& b = loop.points[i];
            if (a.boxDist(b) != 1)
                report("don't expect disconnected loops yet");
            int grow[2] = { EUp, EDown };
            if (a.x == b.x)
            {
                grow[0] = ELeft;
                grow[1] = ERight;
            }

            for (int g = 0; g < 2; g++)
            {
                Pos ga = a + KDir[grow[g]];
                Pos gb = b + KDir[grow[g]];
                if (!placement.isValid(ga) || !placement.isValid(gb))
                    continue;
                if (placement[ga] == n && placement[gb] == n)
                    loop.growth.push_back({ i, {ga, gb} });
            }
        }

        if (loop.start_free && (loop.remaining() % 2 == 1 || loop.growth.empty()))
        {
            const Pos& a = loop.points.front();
            for (int d = 0; d < 4; d++)
            {
                Pos ga = a + KDir[d];
                if (!placement.isValid(ga))
                    continue;
                if (placement[ga] == n)
                {
                    loop.growth.clear();
                    loop.growth.push_back({ 0, {ga} });
                }
            }
        }
    }

    bool test_loop_growth_trouble(Loop loop)
    {
        if (loop.growth.size() <= 1)
            return true;
        vector<Pos> temp_added;
        bool trouble = false;
        while (!trouble && loop.remaining() > 2)
        {
            size_t next_allowed = 0;
            for (const Growth& growth : loop.growth)
            {
                if (loop.remaining() <= 2)
                    break;
                if (growth.where < next_allowed)
                    continue;
                if (growth_ok(growth))
                {
                    add_growth(loop, growth);
                    for (const Pos& p : growth.add)
                        temp_added.push_back(p);
                    for (Growth& later : loop.growth)
                    {
                        if (later.where <= growth.where)
                            continue;
                        later.where += growth.add.size();
                    }
                    next_allowed = growth.where + 1;
                }
            }
            loop.growth.clear();
            add_loop_growth_points(loop);
            if (loop.growth.size() <= 1)
                trouble = true;
        }
        for (const Pos& p : temp_added)
            placement[p] = n;

        return trouble;
    }

    void analyse_loops()
    {
        for (Loop& loop : loops)
        {
            loop.growth.clear();
            loop.single_growth = false;
            if (loop.complete())
                continue;
            loop.pressure.init(N,N,0);
         
            add_loop_growth_points(loop);
            loop.single_growth = test_loop_growth_trouble(loop);

            // pressure map
            double remaining = loop.remaining();
            vector<Pos> now, next;
            if (loop.visited.width == 0)
                loop.visited.init(N, N, 0);
            static int visit_id = 0;
            visit_id++;
            for (const Growth& growth : loop.growth)
            {
                for (const Pos& p : growth.add)
                {
                    if (loop.visited[p] != visit_id)
                        now.push_back(p);
                    loop.visited[p] = visit_id;
                    loop.pressure[p] = remaining;
                }
            }
            remaining -= double(now.size()) * 0.5;
            double pf = 1.0;
            while (remaining > 0 && !now.empty())
            {
                next.clear();
                for (const Pos& p : now)
                {
                    for (int d = 0; d < 4; d++)
                    {
                        Pos q = p + KDir[d];
                        if (!loop.visited.isValid(q))
                            continue;
                        if (loop.visited[q] == visit_id)
                            continue;
                        if (placement[q] != n)
                            continue;
                        next.push_back(q);
                        loop.visited[q] = visit_id;
                        loop.pressure[q] = remaining*pf;
                    }
                }
                now.swap(next);
                remaining -= double(now.size()) * 0.1;
                pf *= 0.9;
            }
        }
    }
    
    void add_growth(Loop& loop, const Growth& growth)
    {
        for (const Pos& p : growth.add)
            placement[p] = loop.id;
        loop.points.insert(loop.points.begin() + growth.where, growth.add.begin(), growth.add.end());
    }

    bool growth_ok(const Growth& growth)
    {
        for (const Pos& p : growth.add)
            if (placement[p] != n)
                return false;
        return true;
    }

    bool grow_single_options()
    {
        bool single_found = false;
        for (Loop& loop : loops)
        {
            if (loop.complete())
                continue;
            if (loop.growth.empty())
                throw PlacementFail(loop.id);
            if (!loop.single_growth)
                continue;
            single_found = true;
            const Growth& growth = loop.growth[0];
            if (!growth_ok(growth))
                throw PlacementFail(loop.id);
            //cout << "single grow " << loop.id << endl;
            add_growth(loop, growth);
        }
        return single_found;
    }
    
    enum CheckResult { Ok, Stop, Undo };

    void grow_highest_pressure()
    {
        Grid<double> pressure(N, N, 0);
        for (Loop& loop : loops)
        {
            if (loop.complete())
                continue;
            FOR_GRID(p, pressure)
                pressure[p] += loop.pressure[p];
        }
        //cout << pressure;

        // todo find highest pressure loop?
        size_t max_r = 0;
        size_t l = 0;
        for (size_t i=0; i<loops.size(); i++)
        {
            const Loop& loop = loops[i];
            size_t r = loop.remaining();
            if (r > max_r)
            {
                max_r = r;
                l = i;
            }
        }

//        Loop& loop = loops[l];
        bool something_added = false;
        for (Loop& loop : loops)
        {
            if (loop.complete())
                continue;
            bool any_found = false;
            Growth best_growth;
            double lowest_pressure = 1e99;
            for (const Growth& growth : loop.growth)
            {
                bool ok = true;
                double gp = 0;
                for (const Pos& p : growth.add)
                {
                    if (placement[p] != n)
                        ok = false;
                    gp += pressure[p];
                }
                if (!ok)
                    continue;
                any_found = true;
                if (gp < lowest_pressure)
                {
                    best_growth = growth;
                    lowest_pressure = gp;
                }
            }
            if (lowest_pressure < 1e99)
            {
                //cout << "pressure grow " << loop.id << " " << best_growth.add[0] << endl;
                for (const Pos& p : best_growth.add)
                    placement[p] = loop.id;
                CheckResult check = check_growth(loop.id);
                if (check == Undo)
                {
                    for (const Pos& p : best_growth.add)
                        placement[p] = n;
                }
                else
                {
                    something_added = true;
                    loop.points.insert(loop.points.begin() + best_growth.where, best_growth.add.begin(), best_growth.add.end());
                    if (check == Stop)
                        break;
                }
            }
            if (!any_found)
                throw PlacementFail(loop.id);
        }
        if (!something_added)
            throw PlacementFail(-1);
    }

    CheckResult check_growth(size_t after_loop_id)
    {
        CheckResult result = Ok;
        for (const Loop& loop : loops)
        {
            if (loop.id <= after_loop_id)
                continue;
            if (loop.complete())
                continue;
            size_t g_count = 0;
            for (const Growth& growth : loop.growth)
            {
                bool ok = true;
                for (const Pos& p : growth.add)
                {
                    if (placement[p] != n)
                        ok = false;
                }
                if (ok)
                    g_count++;
            }
            if (g_count == 0)
                return Undo;
            if (g_count == 1)
                result = Stop;
        }
        return result;
    }

    bool all_loops_done()
    {
        for (Loop& loop : loops)
            if (!loop.complete())
                return false;
        return true;
    }
    
    void place_loop_numbers()
    {
        for (Loop& loop : loops)
        {
            size_t i = loop.first;
            for (const Pos& p : loop.points)
                placement[p] = i++;
        }
    }

    int best_score = 0;
    vector<char> best_path;

    void score_and_record()
    {
        Pos p(N / 2, N / 2);
        size_t i = placement[p];
        vector<char> path;
        path.reserve(i);
        while (i > 0)
        {
            i--;
            bool found = false;
            for (TDir d : {EUp, ERight, EDown, ELeft})
            {
                Pos q = p + KDir[d];
                if (placement[q] == i)
                {
                    found = true;
                    p = q;
                    path.push_back("URDL"[d]);
                }
            }
            if (!found)
                report("path fail");
        }

        int score = problem.score_chars(path);
        if (score > best_score)
        {
            best_score = score;
            best_path = path;
        }
    }
};

vector<char> run(const Problem& p)
{
    SnakeCharmer prog;
    return prog.findSolution(p.N, p.V, p.snake);
}

void tester()
{
    Problem p;
    cin >> p;
    vector<char> ret = run(p);
    cout << ret.size() << endl;
    for (int i = 0; i < (int)ret.size(); ++i)
            cout << ret[i] << endl;
    cout.flush();
}

void capture(const string& save_name)
{
    Problem p;
    cin >> p;
    ofstream out(save_name);
    out << p;
}

void gen(const string& s)
{
    stringstream strm(s);
    int seed;
    strm >> seed;
    Problem p;
    p.generate(seed);
    cout << p;
}

void test(int argc, char** argv)
{
    Problem p = { 7, 3, "3444433322344332432434424422324243233232232424424" };
    //string path = "DLULDDLUUURULURRRDLDRRRULURRDDDDDDLLLLURRUURDD";
    //vector<char> chars(path.begin(), path.end());
    vector<char> chars = run(p);
    cout << p.score_chars(chars) << endl;
}

void runs(int argc, char** argv)
{
    int i = 1;
    if (argc >= 1)
    {
        stringstream strm(argv[0]);
        strm >> i;
    }
    Problem p;
    p.generate(i);

    int v = p.V + 1;
    int n = p.N - 1;
    int c = 4;
    if (argc >= 2)
    {
        stringstream strm(argv[1]);
        strm >> v;
    }
    if (argc >= 3)
    {
        stringstream strm(argv[2]);
        strm >> n;
    }
    if (argc >= 4)
    {
        stringstream strm(argv[3]);
        strm >> c;
    }

    SnakeCharmer prog;
    prog.init(p.N, p.V, p.snake);
    cout << p;
    prog.find_runs(v);
    for (auto r : prog.runs)
        if (r.length > 1) cout << r.start << " " << r.length << endl;
    auto m = prog.best_meeting(n, v, p.N*p.N);
    cout << m.size() << endl;
    for (auto p : m)
        cout << p.first << " " << p.second << endl;
    try
    {
        prog.initial_placement(m, c, v);
        prog.connect_snake();
        cout << "success" << endl;
    }
    catch (PlacementFail f)
    {
        cout << "fail " << f.id << endl;
    }
    Grid<char> dg(p.N, p.N, ' ');
    string ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    FOR_GRID(q, dg)
    {
        size_t x = prog.placement[q];
        if (x == p.N * p.N) continue;
        if (x < 50) dg[q] = ids[x];
        else dg[q] = '#';
    }
    cout << dg;
    for (auto loop : prog.loops)
    {
        cout << ids[loop.id] << loop.id << " "  << loop.first << " " << loop.last << " " << (loop.complete() ? "Y" : "n") << " " << loop.start_free;
        for (Pos p : loop.points)
            cout << " " << p;
        cout << endl;
    }
}

void eval(int argc, char** argv)
{
    int i = 1;
    if (argc)
    {
        stringstream strm(argv[0]);
        strm >> i;
    }
    Problem p;
    p.generate(i);
    cout << p;
    SnakeCharmer prog;
    vector<char> s = prog.findSolution(p.N, p.V, p.snake);
    int score = p.score_chars(s, true);
    cout << "bestV=" << prog.bestV << " bestN=" << prog.bestN << " bestC=" << prog.bestC << endl;
    cout << score << endl;
}

/*
83834
90383
215430
227813
109708
3680
44048
298948
34142
65833
119361
67056
66426
336035
76565
111237
28279
866
44360
155388
*/

void score(int argc, char** argv)
{
    int i = 1, j=1;
    if (argc)
    {
        stringstream strm(argv[0]);
        strm >> i;
        j = i;
    }
    if (argc>1)
    {
        stringstream strm(argv[1]);
        strm >> j;
        if (i > j)
            swap(i, j);
    }
    for (; i <= j; i++)
    {
        Problem p;
        p.generate(i);
        SnakeCharmer prog;
        milliseconds(true);
        vector<char> s = prog.findSolution(p.N, p.V, p.snake);
        int score = p.score_chars(s, false);
        cout << score << endl;
    }
}

int main(int argc, char** argv) 
{
    string arg="tester";
    if (argc>1)
        arg = argv[1];
    
    if (arg=="tester")
        tester();
    else if (arg == "capture" && argc>2)
        capture(argv[2]);
    else if (arg == "gen" && argc > 2)
        gen(argv[2]);
    else if (arg == "test")
        test(argc - 2, argv + 2);
    else if (arg == "runs")
        runs(argc - 2, argv + 2);
    else if (arg == "eval")
        eval(argc - 2, argv + 2);
    else if (arg == "score")
        score(argc - 2, argv + 2);
    else
        cout << "what?" << endl;

    return 0;
}
