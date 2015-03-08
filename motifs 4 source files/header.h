#include <algorithm> 
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
typedef map<string, unsigned int> MotPermFromFile;
typedef map<unsigned int,unsigned int> MotArr;


void Gen_File_Motives();
void Gen_List_of_4_Sus(vector< vector<unsigned int> > &);
void CountMotives_4(vector<unsigned int> &);
