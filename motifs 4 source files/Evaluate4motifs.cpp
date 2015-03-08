#include "stdafx.h"
#include "Evaluate4motifs.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
const mxArray*	pInputGraph	=	NULL; 

unsigned int	NumberOfLinksInGraph	=	0;

//---------------------GLOBALS---------------------------
DirectNetwork directNet;
vector<pair< set<unsigned int>, set<unsigned int> > > directNet1;
set<pair<int,unsigned int>  > Order;
unsigned long array_4[4096];
unsigned long arrMot[4096];
//MotArr MotivesArray;
unsigned int MotivesArray[32000]={0};
//-------------------------------------------------------

void ErrorMsg(char* srt) {
	mxArray *ErrStr = mxCreateString(srt);	
	mxArray *rhs[1] = {ErrStr};
	mxArray *lhs[1];
	mexCallMATLAB(0, lhs, 1, rhs, "disp");
	//mexCallMATLAB(0, lhs, 1, rhs, "pause");
}

void ErrorMsg(char* srt, int i) {
	char A[100];
	char B[100];
	strcpy(A,srt);
	_itoa(i,B,10);
	strcat(A,B);
	mxArray *ErrStr = mxCreateString(A);	
	mxArray *rhs[1] = {ErrStr};
	mxArray *lhs[1];
	mexCallMATLAB(0, lhs, 1, rhs, "disp");
	//mexCallMATLAB(0, lhs, 1, rhs, "pause");
}

void ErrorQuit() {
	char A[100];
	char B[100];
	strcpy(A,"asd");
	_itoa(2,B,10);
	strcat(A,B);
	mxArray *ErrStr = mxCreateString(A);	
	mxArray *rhs[1] = {ErrStr};
	mxArray *lhs[1];
	//mexCallMATLAB(0, lhs, 1, rhs, "disp");
	mexCallMATLAB(0, lhs, 1, rhs, "pause");
}

void	SetInputParameters(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	pInputGraph	=	NULL;

	if (nrhs != 1) { mexErrMsgTxt(" one input arguments required."); }
	if (nlhs != 1){	mexErrMsgTxt("one output arguments allowed.");	}
	pInputGraph	= prhs[0];	
	{
		mxArray *GraphType = mxCreateString("Graph");	
		mxArray *ErrorMessage = mxCreateString("The input must be of the type \"Graph\". Please use ObjectCreateGraph function");	
		mxArray *rhs[3] = {const_cast<mxArray*>(pInputGraph), GraphType, ErrorMessage};
		mxArray *lhs[1];
		mexCallMATLAB(0, lhs, 3, rhs, "ObjectIsType");
		mxDestroyArray(GraphType);		GraphType		=	NULL;
		mxDestroyArray(ErrorMessage);	ErrorMessage	=	NULL;
	}	
	NumberOfLinksInGraph = mxGetM(mxGetField(pInputGraph,0,"Data"));
} 
//-----------------------------------------------------------------------------------------------
void	PrepareToRun(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	const mxArray* const pData = mxGetField(pInputGraph,0,"Data"); 		
	const double* const  DataPtr = mxGetPr(pData); 
	for (unsigned int i = 0; i < NumberOfLinksInGraph; ++i) {
		if (	(unsigned int)floor(*(DataPtr+i+NumberOfLinksInGraph) + 0.5)!= (unsigned int)floor( *(DataPtr+i) + 0.5) ){
			directNet[(unsigned int)floor(*(DataPtr+i+NumberOfLinksInGraph) + 0.5)].first.insert( (unsigned int)floor( *(DataPtr+i) + 0.5));
			directNet[(unsigned int)floor(*(DataPtr+i) + 0.5)].second.insert((unsigned int)floor( *(DataPtr+i+NumberOfLinksInGraph) + 0.5) );
		}
	}
	DirectNetwork::reverse_iterator R=directNet.rbegin();
	char xx[100];
	ErrorMsg(_itoa(R->first,xx,10));
	pair< set<unsigned int>, set<unsigned int> > PTemp;
	for(int i=0;i<=R->first;i++)
		directNet1.push_back(PTemp);
	for(DirectNetwork::iterator D=directNet.begin();D!=directNet.end();++D)
		directNet1[D->first]=D->second;
	ErrorMsg("PrepareToRun Finished");
}
//------------------------------------------------------------------------------------------------------------------ 

vector<vector< unsigned int> > nodesMotifs;

void	PerformCalculations(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{ 

	// Bit strings - matrix of 4x4 excluding the diagonal - 12 bits
	for(int i=0;i<4096;++i)
		array_4[i]=0;

	// Creating the undirected network
	MergedNetwork mergNet;
	for (DirectNetwork::iterator it = directNet.begin(); it != directNet.end(); ++it)
	{
		// inputs
		mergNet[it->first] = it->second.first;
		// outputs one by one
		for (set<unsigned int>::const_iterator it1 = it->second.second.begin(); it1 != it->second.second.end(); ++it1)
			mergNet[it->first].insert(*it1);
	}

	// Now should convert to ordered vector instead of sets (tree)

	// Something like reordering by output degree - not useful
	/*int sum;
	for (MergedNetwork::iterator it = mergNet.begin(); it != mergNet.end(); ++it)
	{
	sum = it->second.size();
	for (set<unsigned int>::iterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
	sum = sum + mergNet[*it1].size();
	Order.insert(std::pair<int,unsigned int>(sum,it->first));
	}*/

	directNet.clear();
	MergedNetwork::reverse_iterator R=mergNet.rbegin();
	// Vector implmentation
	vector< vector<unsigned int> > mergNet1;//(R->first);
	
	// Allocation
	vector<unsigned int> PTemp;
	for(int i=0;i<=R->first;i++)
		mergNet1.push_back(PTemp);
	
	for(MergedNetwork::iterator D=mergNet.begin();D!=mergNet.end();++D) {
		for(set<unsigned int>::iterator S=D->second.begin();S!=D->second.end();++S)
			mergNet1[D->first].push_back(*S);
	}

	mergNet.clear();
	//MotPermFromFile fileMotives;
	Gen_File_Motives();	

	for(int i=1;i<32000;i++)
		MotivesArray[i]=0;
	
	Gen_List_of_4_Sus(mergNet1);

	//Update MotivesArray
	//for(int i=0;i<4096;++i)
	//	MotivesArray[arrMot[i]]+=array_4[i];

	//Free resourses
	mergNet1.clear();
}

MotArr	mapMotives;

void Gen_File_Motives()
{
	ifstream in_file;
	string key;
	unsigned int mot_id, index;

	// ------- Read from file --------------
	in_file.open("4_nodes_data.txt", ios::in);

	while (!in_file.eof())
	{
		in_file >> key >> mot_id >> index;
		arrMot[index] = mot_id;
		mapMotives[mot_id] = 0;
	}


	in_file.close();
	int j = 0;
	for (MotArr::iterator I = mapMotives.begin(); I != mapMotives.end(); ++I){
		I->second = j;
		j++;
	}
	
}

void Gen_List_of_4_Sus(vector< vector<unsigned int> > &merg)
{
	//ErrorMsg("Gen_List_of_4_Sus Start");
	//clock_t start, finish;
	vector<unsigned int> four;
	vector<unsigned int> dummy(200,0);
	vector<unsigned int> StartingPlaces(merg.size(),0);
	for (int i = 0; i < merg.size(); i++){
		StartingPlaces[i] = 0;
		nodesMotifs.push_back(dummy);
	}
	unsigned int first,second,third, forth;
	unsigned int secondIdx, thirdIdx, forthIdx;
//	start=clock();
	int hy;
	//for (set<pair<int,unsigned int> >::reverse_iterator it = Order.rbegin(); it!= Order.rend(); )
	//for (vector<vector<unsigned int> >::iterator it = merg.begin(); it!= merg.end(); 
	ErrorMsg("Gen_List_of_4_Sus 2");

	ErrorMsg("MergSize",merg.size());
	for(first=0;first<merg.size();first++)
	{
		//ErrorMsg("First",first);
		if(merg[first].empty())
			continue;
		//ErrorMsg("FirstOK");
		
		four.push_back(first);//-1-
		for(secondIdx=StartingPlaces[first];secondIdx<merg[first].size();secondIdx++)
		{
			second=merg[first][secondIdx];
			//ErrorMsg("second ",second);
			// Erasing
			StartingPlaces[second]++;

			four.push_back(second);//-2-
			// third is son of first
			for(thirdIdx=secondIdx+1;thirdIdx<merg[first].size();thirdIdx++) {
				third=merg[first][thirdIdx];
				four.push_back(third);//-3-
				// forth is son of first
				for(forthIdx=thirdIdx+1;forthIdx<merg[first].size();forthIdx++) {
					// Father and 3 sons
					four.push_back(merg[first][forthIdx]);
					CountMotives_4(four);
					four.pop_back(); // forth
				}

				// Now Motif like Grandfather, father, son, left uncle
				// i.e. Grandfather, 2 fathers, son of right father
				for(forthIdx=StartingPlaces[third]+1;forthIdx<merg[third].size();forthIdx++) {
					forth=merg[third][forthIdx];
					bool ForthIsGrandFather = forth==first;
					bool ForthIsUncle = !(find(merg[first].begin(), merg[first].end(), forth)==merg[first].end());
					// Checks that forth is not the mutual son of the two brothers
					bool ForthIsIncest = !(find(merg[second].begin(), merg[second].end(), forth)==merg[second].end());
					bool ForthIsProblem = ForthIsUncle || ForthIsIncest || ForthIsGrandFather;
					if(!ForthIsProblem) {
						four.push_back(forth);
						CountMotives_4(four);
						four.pop_back(); // forth
					}
				}
				if(four.size()!=3) {
					ErrorMsg("Not thirdX"); return;	}
				four.pop_back(); // third
			}
			// Third is son of second
			for(thirdIdx=StartingPlaces[second];thirdIdx<merg[second].size();thirdIdx++) {
				third=merg[second][thirdIdx];
				bool ThirdIsUncle = !(find(merg[first].begin(), merg[first].end(), third)==merg[first].end());
				bool ThirdAlreadyComputed = third<first;
				bool ThirdProblem = ThirdIsUncle || ThirdAlreadyComputed;
				//ErrorMsg("Third",third);
				//ErrorMsg("ThirdIsUncle",ThirdIsUncle);
				//ThirdIsUncle=false;
				if(  !ThirdProblem   ) {
					four.push_back(third);
					// Now Motif like Grandfather, father, 2 sons
					for(forthIdx=thirdIdx+1;forthIdx<merg[second].size();forthIdx++) {
						forth=merg[second][forthIdx];
						bool ForthIsGrandFatherOrComputed = forth<=first;
						bool ForthIsUncle = !(find(merg[first].begin(), merg[first].end(), forth)==merg[first].end());
						bool ForthProblem = ForthIsUncle || ForthIsGrandFatherOrComputed;
						if( !ForthProblem) 	
						{
							// 1 GrandFather, 2 - Father, and now looking for 2 grand sons which are not uncles
							four.push_back(forth);
							CountMotives_4(four);
							four.pop_back();			
						}
					}
					// Now Motif like Grandfather, father, son, grandson
					for(forthIdx=0;forthIdx<merg[third].size();forthIdx++) {
						forth=merg[third][forthIdx];
						bool ForthIsGrandFatherOrComputed = forth<=first;
						bool ForthIsGreatUncle = !(find(merg[first].begin(), merg[first].end(), forth)==merg[first].end());
						bool ForthIsUncle = !(find(merg[second].begin(), merg[second].end(), forth)==merg[second].end());
						bool ForthProblem = ForthIsGreatUncle || ForthIsUncle || ForthIsGrandFatherOrComputed;
						if( !ForthProblem    )
						{
							four.push_back(forth);
							CountMotives_4(four);
							four.pop_back();
						}
					}

					// Now Motif like Grandfather, father, son, right uncle
					//for (set<unsigned int>::const_iterator iter2 = iter1; iter2 != merg[second].end(); ++iter2)//check all the other grand-children
					for(forthIdx=secondIdx+1;forthIdx<merg[first].size();forthIdx++) {
						//ErrorMsg("XX74");
						//ErrorMsg("XX",first);
						//ErrorMsg("XX",second);
						//ErrorMsg("XX",third);
						//ErrorMsg("XX",forth);

						forth=merg[first][forthIdx];
						// Mistake here
						//bool ForthIsGreatUncle = find(merg[first].begin(), merg[first].end(), forth)==merg[first].end();
						// forth is not third
						//if(!ForthIsGreatUncle) {
						if(forth!=third) {
							four.push_back(forth);
							CountMotives_4(four);
							four.pop_back();
						}
					}
				if(four.size()!=3) {
					ErrorMsg("Not third"); return;	}
				four.pop_back(); // third
				}
			}
			if(four.size()!=2) {
					ErrorMsg("Not second"); return;	}
			four.pop_back(); // second
		}
		if(four.size()!=1) {
					ErrorMsg("Not 1"); return;	}
		four.pop_back(); // first
		// Assert four is empty
		if(!four.empty()) {
			ErrorMsg("Not empty XJXJX");	}
	}
	merg.clear();
}
map<vector<unsigned int>,unsigned int> DebugAllFours;
void CountMotives_4(vector<unsigned int> & four)
{
	unsigned int first = four[0];
	unsigned int second = four[1];
	unsigned int third = four[2];
	unsigned int fourth = four[3];
	vector<unsigned int> TVec = four;
	sort(TVec.begin(), TVec.end());
	//if (DebugAllFours.find(four) == DebugAllFours.end()) {
	if (DebugAllFours[TVec] == 0) {
		//DebugAllFours[four]=1;
		DebugAllFours[TVec] = four[0] + four[1] * 100 + four[2] * 10000 + four[3] * 1000000;
	}
	else {
		ErrorMsg("DD", DebugAllFours[TVec]);
		ErrorMsg("DebugAllFours-----");
		ErrorMsg("DebugAllFours", four[0]);
		ErrorMsg("DebugAllFours", four[1]);
		ErrorMsg("DebugAllFours", four[2]);
		ErrorMsg("DebugAllFours", four[3]);
		ErrorQuit();
	}

	int temp;
	//set<unsigned int>::iterator I=four.begin();
	temp = 0;

	//---------------- 1 ----------------------------
	if (directNet1[first].second.find(second) != directNet1[first].second.end())
		temp = temp + 2048;
	if (directNet1[first].second.find(third) != directNet1[first].second.end())
		temp = temp + 1024;
	if (directNet1[first].second.find(fourth) != directNet1[first].second.end())
		temp = temp + 512;
	//---------------- 2 ----------------------------
	if (directNet1[second].second.find(first) != directNet1[second].second.end())
		temp = temp + 256;
	if (directNet1[second].second.find(third) != directNet1[second].second.end())
		temp = temp + 128;
	if (directNet1[second].second.find(fourth) != directNet1[second].second.end())
		temp = temp + 64;
	//---------------- 3 ----------------------------
	if (directNet1[third].second.find(first) != directNet1[third].second.end())
		temp = temp + 32;
	if (directNet1[third].second.find(second) != directNet1[third].second.end())
		temp = temp + 16;
	if (directNet1[third].second.find(fourth) != directNet1[third].second.end())
		temp = temp + 8;
	//---------------- 4 ----------------------------
	if (directNet1[fourth].second.find(first) != directNet1[fourth].second.end())
		temp = temp + 4;
	if (directNet1[fourth].second.find(second) != directNet1[fourth].second.end())
		temp = temp + 2;
	if (directNet1[fourth].second.find(third) != directNet1[fourth].second.end())
		temp = temp + 1;


	array_4[temp]++;
	
	/*
	MotArr MotivesArray1;
	*/
	for (int i = 0;i<4 ; i++){
		nodesMotifs[four[i]][mapMotives[arrMot[temp]]]++;
	}

	
}

//================================================================================================
void	FreeResourses(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	MotArr MotivesArray;
	for(int i=0;i<4096;++i)
		MotivesArray[arrMot[i]]=0;
	//Update MotivesArray
	for(int i=0;i<4096;++i)
		MotivesArray[arrMot[i]]+=array_4[i];

	mxArray* motifs4_length=mxCreateDoubleMatrix(MotivesArray.size(),2,mxREAL);
	double* P = mxGetPr(motifs4_length);
	
	for (std::map<unsigned int,unsigned int>::const_iterator LengthIt = MotivesArray.begin();LengthIt !=  MotivesArray.end(); LengthIt++) {
		*P=LengthIt->first;
		*(P+MotivesArray.size())=LengthIt->second;
		++P;
	}
	plhs[0] =  motifs4_length;
	
	ofstream myfile;
	myfile.open("example.txt");
	for (int i = 0; i < nodesMotifs.size(); ++i){
		for (int j = 0; j < MotivesArray.size() ; ++j)
			myfile << nodesMotifs[i][j] << "\t";
		myfile << endl;
		}
	myfile.close();
	
	
	
	//delete P;
	
	directNet1.clear();
	MotivesArray.clear();
	nodesMotifs.clear();
	//Order.clear();
}



