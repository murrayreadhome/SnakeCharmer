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

using namespace std;

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
    for (size_t i=0; i<N*N; i++)
      snake[i] = '0' + uniform_int_distribution<int>(2,V+1)(re);
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

class SnakeCharmer {
public:
  vector<char> findSolution(int N, int V, string snake)
  {           
    char names[] = {'L','D','R','U'};  

    int n=N*N-1;
    vector<char> moves(n);
    
    for (int i=0,dir=0,L=1; i<n; )
    {
      for (int k=0; k<L && i<n; k++,i++) moves[i]=names[dir];
      dir=(dir+1)%4;
      for (int k=0; k<L && i<n; k++,i++) moves[i]=names[dir];
      dir=(dir+1)%4;
      L++;
    }
    
    return moves;
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

int main(int argc, char** argv) 
{
  string arg="tester";
  if (argc>1)
    arg = argv[1];
  
  if (arg=="tester")
    tester();
  else if (arg == "capture" && argc>2)
    capture(argv[2]);
  else if (arg == "gen" && argc>2)
    gen(argv[2]);

  return 0;
}
