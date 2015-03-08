#include <vector> 
#include <iostream> 
#include <map>
#include <list>
#include <fstream>
#include <set>
#include <string>

using namespace std;


typedef map<unsigned int, pair< set<unsigned int>, set<unsigned int> > > DirectNetwork;
typedef map<unsigned int, set<unsigned int> > MergedNetwork;
typedef set<vector<unsigned int> > Suspected;
typedef map<int,int> MotArr;

void CountMotives(vector<unsigned int> &);


