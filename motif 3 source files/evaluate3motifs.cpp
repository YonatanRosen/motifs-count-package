#include "stdafx.h"
#include "Evaluate3motifs.h"
#include<math.h>

const mxArray*	pInputGraph	=	NULL; 

unsigned int	NumberOfLinksInGraph	=	0;
MotArr	mapMotives;

void ErrorMsg(char* srt) {
	mxArray *ErrStr = mxCreateString(srt);
	mxArray *rhs[1] = { ErrStr };
	mxArray *lhs[1];
	mexCallMATLAB(0, lhs, 1, rhs, "disp");
	//mexCallMATLAB(0, lhs, 1, rhs, "pause");
}

void ErrorMsg(char* srt, int i) {
	char A[100];
	char B[100];
	strcpy(A, srt);
	_itoa(i, B, 10);
	strcat(A, B);
	mxArray *ErrStr = mxCreateString(A);
	mxArray *rhs[1] = { ErrStr };
	mxArray *lhs[1];
	mexCallMATLAB(0, lhs, 1, rhs, "disp");
	//mexCallMATLAB(0, lhs, 1, rhs, "pause");
}

void ErrorQuit() {
	char A[100];
	char B[100];
	strcpy(A, "asd");
	_itoa(2, B, 10);
	strcat(A, B);
	mxArray *ErrStr = mxCreateString(A);
	mxArray *rhs[1] = { ErrStr };
	mxArray *lhs[1];
	//mexCallMATLAB(0, lhs, 1, rhs, "disp");
	mexCallMATLAB(0, lhs, 1, rhs, "pause");
}
//---------------------GLOBALS---------------------------
DirectNetwork directNet;
Suspected list_of_sus;
int array_3[64];
int arrMot[64];
MergedNetwork mergNet;
MotArr MotivesArray;
ofstream out_file;
vector<vector< unsigned int> > nodesMotifs;

//------------------------------------------------------------------------------------------------------------------ 

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
		
}
//------------------------------------------------------------------------------------------------------------------ 
void	PerformCalculations(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{ 
	for(int i=0;i<64;++i)
		array_3[i]=0;

// #####  MergeNetLists  #####
	for (DirectNetwork::iterator itd = directNet.begin(); itd != directNet.end(); ++itd)
	{
		mergNet[itd->first] = itd->second.first;

		for (set<unsigned int>::const_iterator it1d = itd->second.second.begin(); it1d != itd->second.second.end(); ++it1d)
			mergNet[itd->first].insert(*it1d);
	}
// #####  End of MergeNetLists #####
	vector<unsigned int> dummy(14, 0);
	vector<unsigned int> StartingPlaces(directNet.size() + 1, 0);
	for (int i = 0; i <= directNet.size(); i++){
		StartingPlaces[i] = 0;
		nodesMotifs.push_back(dummy);
	}

//  ##### Gen_File_Motives #####
	ifstream in_file;
	string key;
	unsigned int mot_id, index;

	in_file.open("3_nodes_data.txt",  ios :: in);

	while ( !in_file.eof() )
	{
		in_file>>key>>mot_id>>index;
		arrMot[index] = mot_id;//int arrMot[64];
		mapMotives[mot_id] = 0;
	}
	int j = 0;
	for (MotArr::iterator I = mapMotives.begin(); I != mapMotives.end(); ++I){
		I->second = j;
		j++;
	}
	in_file.close();
//  ##### End of Gen_File_Motives #####

	
//  #####  InitMotivesArray #####
	

	MotivesArray[6] = 0;
	MotivesArray[12] = 0;
	MotivesArray[14] = 0;
	MotivesArray[36] = 0;
	MotivesArray[38] = 0;
	MotivesArray[46] = 0;
	MotivesArray[74] = 0;
	MotivesArray[78] = 0;
	MotivesArray[98] = 0;
	MotivesArray[102] = 0;
	MotivesArray[108] = 0;
	MotivesArray[110] = 0;
	MotivesArray[238] = 0;
//  #####  End of InitMotivesArray #####

//  #####  Gen_List_of_3_Sus  #####
	vector<unsigned int> threes;
	out_file.open("debugeer.txt" , ios::out);

	for (MergedNetwork::const_iterator it = mergNet.begin(); it != mergNet.end(); )
	{
		threes.push_back(it->first);//first
		for (set<unsigned int>::const_iterator it1 = mergNet[it->first].begin(); it1 != mergNet[it->first].end(); ++it1)
		{
			
				threes.push_back(*it1);//second
			
			set<unsigned int>::const_iterator help = it1;
			help++;
			for (;help != mergNet[it->first].end(); ++help)
			{
			
				threes.push_back(*help);//third
				CountMotives(threes);
				threes.pop_back(); //pop third
			}		
			threes.pop_back();//pop second
			
		}

		for (set<unsigned int>::const_iterator iter1 = mergNet[it->first].begin(); iter1 != mergNet[it->first].end(); ++iter1)
		{
			threes.push_back(*iter1);//second
			mergNet[*iter1].erase(mergNet[*iter1].begin());
			for (set<unsigned int>::const_iterator iter2 = mergNet[*iter1].begin(); iter2 != mergNet[*iter1].end(); ++iter2)
			{
				if ( (*iter2 != it->first ) && ( mergNet[it->first].find(*iter2)==mergNet[it->first].end())){
					threes.push_back(*iter2);
					CountMotives(threes);
					threes.pop_back(); //pop third
				}
			}
			threes.pop_back();
		}
	threes.pop_back();

		++it;
		if (it != mergNet.end())
			mergNet.erase(mergNet.begin());
	}
//  #####  End of Gen_List_of_3_Sus  #####

    out_file.close();
//  #####  UpdateMotivesArray  #####
	for(int i=0;i<64;++i)
		MotivesArray[arrMot[i]]+=array_3[i];
//  #####  End of UpdateMotivesArray  #####
}

void CountMotives(vector<unsigned int> &threes)
{
	int temp;
	set<unsigned int>::iterator it_find;
	
//	for (Suspected::iterator it = list_of_sus.begin(); it != list_of_sus.end(); ++it)
//	{
		unsigned first = threes[0];
		unsigned int second = threes[1];
		unsigned int third = threes[2];
		temp=0;
		
		
		//---------------- 1 ----------------------------
		it_find = directNet[first].second.find(second);
		if (it_find != directNet[first].second.end())
			temp = temp+32;
		it_find = directNet[first].second.find(third);
		if (it_find != directNet[first].second.end())
			temp = temp+16;
		//---------------- 2 ----------------------------
		it_find = directNet[second].second.find(first);
		if (it_find != directNet[second].second.end())
			temp = temp+8;
		it_find = directNet[second].second.find(third);
		if (it_find != directNet[second].second.end())
			temp = temp+4;
		//---------------- 3 ----------------------------
		it_find = directNet[third].second.find(first);
		if (it_find != directNet[third].second.end())
			temp = temp+2;
		it_find = directNet[third].second.find(second);
		if (it_find != directNet[third].second.end())
			temp = temp+1;
		//------- Finished analyzing the sus. three --------------
	//if (temp==25 || temp==38)
	//	out_file<<first<<" "<<second<<" "<<third<<endl;
		array_3[temp]++;//Updating the 'int array_3[64];'

		for (int i = 0; i<3; i++){
			nodesMotifs[threes[i]][mapMotives[arrMot[temp]]]++;
		}
	//}
}
//================================================================================================
void	FreeResourses(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{  

	mxArray* motifs3_length=mxCreateDoubleMatrix((int)MotivesArray.size(),2,mxREAL);
	double* P = mxGetPr(motifs3_length);


	
	for (std::map<int,int>::const_iterator LengthIt = MotivesArray.begin();LengthIt !=  MotivesArray.end(); LengthIt++) {
		*P=LengthIt->first;
		*(P+MotivesArray.size())=LengthIt->second;
		++P;
	}
	ofstream myfile;
	myfile.open("example3.txt");
	for (int i = 0; i < nodesMotifs.size(); ++i){
		for (int j = 0; j < MotivesArray.size(); ++j)
			myfile << nodesMotifs[i][j] << "\t";
		myfile << endl;
	}
	myfile.close();
	/*for (std::map<unsigned int,double >::const_iterator LengthCircleIt = Circles.begin();LengthCircleIt !=  Circles.end(); LengthCircleIt++) {
	numberOfCircles=(int)LengthCircleIt->second;
		*P=LengthCircleIt->second;
		++P;
	}*/
	plhs[0] =  motifs3_length;
	
	directNet.clear();
	list_of_sus.clear();
	mergNet.clear();
	MotivesArray.clear();

}



