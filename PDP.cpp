#include <cstring>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stack>
#include <stdio.h>
#include <time.h> 
using namespace std;
struct Snapshot{
	vector<int> L1; 
    vector<int> X1; 
	vector< vector<int> > solutions;
}; 
bool BBd_=false, BBb_=false, BBb2_=false;
int problem_size;

//Max element in the input is unique
int width;
vector< vector<int> > solutions;
int solutionsSize;
char* fileName_time= new char[100];
char* fileName_output = new char[100];

ofstream writefile_time,writefile_output;



bool contains(vector<int> L, int key){
	for(int i=0;i<L.size();i++)
		if(L.at(i)==key)
			return true;
	return false;
}
int findMaxPos(vector<int> L) {

		int maxDigest = 0;
        int maxpos = 0;
		for (int i = 0; i < L.size(); i++) {
			if (L.at(i)>= maxDigest) {
				maxDigest = L.at(i);
                maxpos = i;
			}
		}
		return maxpos;

    }
bool contains(vector< vector<int> > L, vector<int> key){
	for(int i=0;i<L.size();i++)
		if(L.at(i)==key)
			return true;
	return false;
}
void copy_o(vector<int> d,vector<int> &d1){
           for (int i = 0; i < d.size(); i++) {
			   d1.push_back(/*.insert(*/d.at(i));
			
            }
	//d1 = d;
            
}

vector<int> convert_arry_to_list(int tt[],int n){
        vector<int> L;
        for (int e=0; e < n;e++)
			L.push_back(tt[e]);
        return L;
}



vector<int> getDelta(int partialDigestDistance,vector<int> X) {

		vector<int> resultList; // = new vector<int>();

		for (int ii =0;ii<X.size();ii++){// : X) {
			int rs = X[ii];
			int distance = abs(rs - partialDigestDistance);
			resultList.push_back(distance);
		}
		return resultList;

	}
int find_pos(vector<int> ff, int key){
		int index = -1;
		for (int i=0;i<ff.size();i++){
			if(ff.at(i)== key){
				index = i;
				break;
			}

		}
		return index;
}






























bool isDistanceSubsetOf(vector<int> distanceList,vector<int> L) {

		bool result = true;
		int uu;
		int index;
		for (int i=0;i<distanceList.size();i++){
			
			uu = distanceList.at(i);
			if (!contains(L,uu)) {
				result = false;
				break;
			}
			else
			{
				index = find_pos(L,uu);
				L.erase(L.begin()+index);
			}
		}
		return result;

}
vector<int> removePartialDigests(vector<int> delta, vector<int> &L){
            vector<int> result; // = new ArrayList<Integer>();
			int uu,index;
			for (int i = 0; i < delta.size(); i++) {
                        uu = delta.at(i);//.get(i);
                        index= find_pos(L,uu);
                     	L.erase(L.begin()+index);//.remove(index);
						//L.erase(remove(index,index,uu),L.end());
                        result.push_back(uu);//add(uu);
			}
			return result;
            
}


void print_vector(vector<int> dd){
	cout<<"{";
		for(int j=0;j<dd.size();j++)
		{
			cout<<dd.at(j)<<",";
		}
		cout<<"},\t";
		
}

void print_vector_in_file(vector<int> dd,char* fileName_output)
	{
	writefile_output.open(fileName_output,ios::app);
	writefile_output<<"{";
		for(int j=0;j<dd.size();j++)
		{
			writefile_output<<dd.at(j)<<",";
		}
		writefile_output<<"}\n";
		writefile_output.close();
		
}





void print_solutions_in_file (vector< vector<int> > solutions, char* fileName_output)
{
	writefile_output.open(fileName_output,ios::app);
	writefile_output<<"\nsolutions = {";
	
	writefile_output.close();
	for (int i =0;i<solutions.size();i++)
	{
		print_vector_in_file(solutions.at(i),fileName_output);
	}
	writefile_output.open(fileName_output,ios::app);
	cout<<"}\n";
	writefile_output.close();
}



void print_solutions (vector< vector<int> > solutions){
	cout<<"solutions = {";
	for (int i =0;i<solutions.size();i++)
	{
		print_vector(solutions.at(i));
	}
	cout<<"}\n";
}

void ReadSet(char *szFile, vector <int> &x, int &N)
{
	ifstream myfile;
	myfile.open(szFile, std::ios_base::in);
	if(myfile)
	{
		int a;
		while (myfile >> a)
		{
			//printf("%f", a);
			x.push_back(a);
			N++;
		}
	}
	else
	{
		cout << "An exception occurred. I/O Exception . The file " <<szFile<<" not exist"   << '\n';
		exit(1);
	}

}



bool compare_two_vectors(vector<int> u,vector<int> v)
{
	if(u.size()!=v.size())
		return false;
	sort(u.begin(),u.end());
	sort(v.begin(),v.end());
	for(int i=0;i<u.size();i++)
		if(u[i]!=v[i])
			return false;
	return true;
}
bool compare_two_snapshots(Snapshot e1, Snapshot e2)
{
	//bool resultL=compare_two_vectors(e1.L1,e2.L1);
	//bool resultX=compare_two_vectors(e1.X1,e2.X1);
	if(compare_two_vectors(e1.X1,e2.X1))
	{
		//if(compare_two_vectors(e1.L1,e2.L1))
			return true;
	}
	/*if(resultL&& resultX)
		return true;*/
	return false;

}


bool compare_element_to_vector(vector<Snapshot> PP,Snapshot e)
{
	for(int i=PP.size()-1;i>=0;i--)
		if(compare_two_snapshots(PP[i],e))
			return true;
	return false;
}

vector<Snapshot> divide(vector<Snapshot> snapshot_stack)
{
	Snapshot currentSnapshot;
	vector<int> deltaList_1,deltaList_2;
	bool tst1,tst2;
	int partialDigestDistance,y_pos;
	vector<Snapshot> result;
	result.clear();
	for(int i=0;i<snapshot_stack.size();i++)
	{
			currentSnapshot = snapshot_stack[i];

			if (currentSnapshot.L1.empty()) { // 1 - 3
				sort(currentSnapshot.X1.begin(),currentSnapshot.X1.end());
				
				if(!contains(solutions,currentSnapshot.X1))
					solutions.push_back(currentSnapshot.X1);
				continue;
			}
			y_pos = findMaxPos(currentSnapshot.L1); // 4
			partialDigestDistance = currentSnapshot.L1.at(y_pos);
			deltaList_1 = getDelta(partialDigestDistance, currentSnapshot.X1);
			deltaList_2 = getDelta(width-partialDigestDistance, currentSnapshot.X1);
			tst1 = isDistanceSubsetOf(deltaList_1, currentSnapshot.L1);
			tst2 = isDistanceSubsetOf(deltaList_2, currentSnapshot.L1);
			if(!tst1 && !tst2)
				continue;
			else if (tst1 && tst2)
			{
				
				vector<int> L2; 
				vector<int> X2; 
				copy_o(currentSnapshot.L1,L2);
				copy_o(currentSnapshot.X1,X2);
				currentSnapshot.X1.push_back(partialDigestDistance);
				vector<int> removed = removePartialDigests(deltaList_1, currentSnapshot.L1);
				if(!compare_element_to_vector(result,currentSnapshot))
					result.push_back(currentSnapshot);
				currentSnapshot.L1 = L2;
				currentSnapshot.X1 = X2;
				currentSnapshot.X1.push_back(width-partialDigestDistance);
				vector<int> removed_1 = removePartialDigests(deltaList_2, currentSnapshot.L1);
				if(!compare_element_to_vector(result,currentSnapshot))
					result.push_back(currentSnapshot);
			}
			else if (tst1)
			{
				
				currentSnapshot.X1.push_back(partialDigestDistance);
				vector<int> removed = removePartialDigests(deltaList_1, currentSnapshot.L1);
				if(!compare_element_to_vector(result,currentSnapshot))
					result.push_back(currentSnapshot);
			}
			else if (tst2)
			{
				
				currentSnapshot.X1.push_back(width-partialDigestDistance);
				vector<int> removed_1 = removePartialDigests(deltaList_2, currentSnapshot.L1);
				if(!compare_element_to_vector(result,currentSnapshot))
					result.push_back(currentSnapshot);
			}
	}

	
	return result;
}






void place_BBb(vector<int> L, vector<int> X){

    
	vector<Snapshot> snapshot_stack;
	Snapshot currentSnapshot;
	currentSnapshot.L1 = L;
	currentSnapshot.X1 = X;
	currentSnapshot.solutions = solutions;
	snapshot_stack.push_back(currentSnapshot);

	while(!snapshot_stack.empty())
		snapshot_stack = divide(snapshot_stack);
	
}
int find_minimum(vector<long double> Mem){
	int index = 0;
	long double minimum = Mem[0];
	for(int Alph_0=1;Alph_0<Mem.size();Alph_0++)
	{
		if (Mem[Alph_0]<=minimum)
		{
			minimum = Mem[Alph_0];
			index = Alph_0;
		}
	}
	return index+1;

}
long double choose(const long n, const long r) {
  if(r>n)
    return 0;

  if(r==0||r==n)
    return 1;
  
  long double total=(long double)n/r;
  for(long i=1;i<r;i++)
    total=total*(n-i)/(r-i);
  return total;
}
int find_minimum_memory_mod_13_mod(int n)
{
	vector<long double> Mem;
	for(int i =1;i<=n;i++)
	{
		long double x;
		Mem.push_back(x);
	}
	for(int i=0;i<Mem.size();i++)
	{
		long double sum =0.0;
		for(int j= 1;j<=i+1;j++)
			sum=sum+j;
		Mem[i]= pow((double)2,i)*(i+2+choose(n+1,2)-sum)+pow((double)2,n-i)*(n+2);
	}
	int level_with_minimum_memory = find_minimum(Mem);
	return level_with_minimum_memory;
}

void divide_till_end(vector<Snapshot> snapshot_stack){
	while(!snapshot_stack.empty())
		snapshot_stack = divide(snapshot_stack);
}

void place_BBb2(vector<int> L, vector<int> X){

    
	int limit;
	vector<Snapshot> snapshot_stack;
	Snapshot currentSnapshot;
	currentSnapshot.L1 = L;
	currentSnapshot.X1 = X;
	currentSnapshot.solutions = solutions;
	snapshot_stack.push_back(currentSnapshot);
	int No_of_steps = 0 ;
	limit = find_minimum_memory_mod_13_mod(problem_size);
	cout<<"limit =					"<<limit<<"\n";
	for(int xc=1;xc<=limit;xc++)
	{
		
		if(!snapshot_stack.empty())
			No_of_steps++;
		snapshot_stack = divide(snapshot_stack);
		cout<<"StackSize = "<<snapshot_stack.size()<<"\n";

	}
	cout<<"NoOfSteps = "<<No_of_steps<<"\n";
	/*
	while(!snapshot_stack.empty())
	{
		snapshot_stack = divide(snapshot_stack);
		if(snapshot_stack.size()>=(int)pow(2,limit))
			break;

	}*/
	vector<Snapshot> temp_stack;
	for(int ii=0;ii<snapshot_stack.size();ii++)
	{
		temp_stack.clear();
		temp_stack.push_back(snapshot_stack[ii]);
		divide_till_end(temp_stack);
	}
}



vector<int> MakeHardInstance(int n)
    {
        //Previous implementation
        vector<int> R;
        vector<int> res;
        int i, k, t, n0;
        int m = 1000;  
        int epsilon;
        
        n0 = n / 5;
        R.clear();
		epsilon =rand()%(m/2) +1;
        
        //A2
        for (i = 1; i <= n0; i++)
            {
                t = (i * epsilon);
				if (find_pos(R,t) < 0)
				{ 
					R.push_back(t);
				}
             }
        
        //A4
        for (i = 1; i <= n0; i++)
            {
                t = (2 * n0 + i)*epsilon;
                if (find_pos(R,t) < 0)
                { 
					R.push_back(t);
				}
            }
        
        //A5
        for (i = 1; i <= n0; i++)
            {
                t = n*m - (2 * n0 + i) * epsilon;
                if (find_pos(R,t) < 0)
                { 
					R.push_back(t);
				}
            }
        
        //A1
        for (i = 1; i <= n0; i++)
            {
                t = n*m - i * epsilon;
                if (find_pos(R,t) < 0)
                { 
					R.push_back(t);
				}
            }
        
        //D
        for (i = 1; i <= n0; i++)
          {
            
			  k= rand() %2;
            if (k == 0)
            {
                t = (n0 + i) * epsilon;
            }
            else
            {
                t = n*m - (n0 + i)*epsilon;
            }
            if (find_pos(R,t) < 0)
                { 
					R.push_back(t);
				}
          }
		R.push_back(n*m);
        sort(R.begin(),R.end()); 
        return (R);
}
  
vector<int> MakeRandomInstance(int n, int w)
    {

		vector<int> R;
        int rr;
		//srand (time(NULL));
		while (R.size() < n)
        {
            
			rr = rand()%w +1;
            if (contains(R,rr) == false)
				R.push_back(rr);
        }
        
		sort(R.begin(),R.end());


        
        return (R);
    }
vector<int> ComputeDistance(vector<int> A)
    {
        vector<int> R; //= new vector<int>();
        int k, t, s, m;
        s = 0;
        t = 1;
        m = (A.size() * (A.size() + 1)) / 2;
        for (k = 0; k < m; k++)
        {
            if (k < A.size())
            {
				R.push_back(A.at(k));
                continue;
            }
            if (s + t == A.size())
            {
                s = 0;
                t++;
            }
            
			R.push_back(abs((int)A.at(s+t)-(int)A.at(s)));
            s++;
        }
		sort(R.begin(),R.end());
        return (R);
}


void place_BBd(vector<int> L, vector<int> X){

     
	stack<Snapshot> snapshot_stack;
	Snapshot currentSnapshot;
	currentSnapshot.L1 = L;
	currentSnapshot.X1 = X;
	currentSnapshot.solutions = solutions;
	snapshot_stack.push(currentSnapshot);
	vector<int> deltaList_1,deltaList_2;
	bool tst1,tst2;
	int partialDigestDistance,y_pos;
	while(!snapshot_stack.empty()){
		//cout<<"Stack Size	=	"<<snapshot_stack.size()<<"\n";
		currentSnapshot = snapshot_stack.top();
		snapshot_stack.pop();
		

		if (currentSnapshot.L1.empty()) { // 1 - 3
            sort(currentSnapshot.X1.begin(),currentSnapshot.X1.end());
			#pragma omp critical
			if(!contains(solutions,currentSnapshot.X1))
				solutions.push_back(currentSnapshot.X1);
			continue;
		}
		y_pos = findMaxPos(currentSnapshot.L1); // 4
		partialDigestDistance = currentSnapshot.L1.at(y_pos);
		deltaList_1 = getDelta(partialDigestDistance, currentSnapshot.X1);
		deltaList_2 = getDelta(width-partialDigestDistance, currentSnapshot.X1);
		tst1 = isDistanceSubsetOf(deltaList_1, currentSnapshot.L1);
		tst2 = isDistanceSubsetOf(deltaList_2, currentSnapshot.L1);
		if(!tst1 && !tst2)
			continue;
		else if (tst1 && tst2)
		{
			vector<int> L2; 
			vector<int> X2; 
			copy_o(currentSnapshot.L1,L2);
			copy_o(currentSnapshot.X1,X2);
			
		
		 	currentSnapshot.X1.push_back(partialDigestDistance);
			vector<int> removed = removePartialDigests(deltaList_1, currentSnapshot.L1);
			snapshot_stack.push(currentSnapshot);
			currentSnapshot.L1 = L2;
			currentSnapshot.X1 = X2;
			currentSnapshot.X1.push_back(width-partialDigestDistance);
			vector<int> removed_1 = removePartialDigests(deltaList_2, currentSnapshot.L1);
			snapshot_stack.push(currentSnapshot);
		}
		else if (tst1)
		{
			currentSnapshot.X1.push_back(partialDigestDistance);
			vector<int> removed = removePartialDigests(deltaList_1, currentSnapshot.L1);
			snapshot_stack.push(currentSnapshot);
		}
		else if (tst2)
		{
			currentSnapshot.X1.push_back(width-partialDigestDistance);
			vector<int> removed_1 = removePartialDigests(deltaList_2, currentSnapshot.L1);
			snapshot_stack.push(currentSnapshot);
		}
		
	}
}




void partialDigestProblem(vector<int> L){
		int width_pos = findMaxPos(L); // 1
		clock_t begin,end;
		double start,end_1;
		double elapsed_secs;
		width = L.at(width_pos);      
        cout<<"width = "<<width<<"\n";
		L.erase(L.begin()+width_pos);
        vector<int> X;
		X.push_back(0);
		X.push_back(width);
		vector<int> X1 = X;
		vector<int> L1 = L;
		vector<int> X2 = X;
		vector<int> L2 = L;
		vector< vector<int> > solutions1,solutions2,solutions3,solutions4;
		cout<<"\n\n\n";
		
				
		solutions.clear();
		if(BBd_)
			place_BBd(L2,X2);
		else if (BBb_)
			place_BBb(L2,X2);
		else 
			place_BBb2(L2,X2);

		
		
		
		solutions4 = solutions;
		print_solutions(solutions4);
		solutionsSize = solutions4.size();
		//cout<<fileName_output<<"\n";		
		print_solutions_in_file(solutions4,fileName_output);
		//solutions = solutions4;
		//cout<<"\n\n\n";
		solutions=solutions4;
		
		
}
 




int main(int argc, char* argv[]) {
	int d_selection,a_selection;
	int n,m;
	
	vector <int> rand_solution;
	vector <int> rand_instance;
	vector < vector<int> > solutions;
	double runningTime;
	strcpy(fileName_time,"pdp_time.txt");

	char* Algorithm= new char[100];
if(argc!=2)
	
{
	cout<<"If you want to run the program on specific data; please rerun the programme using the prammeteres like:  \nPDP input_file\n In case of specific data, we will run the BBb2 algorithm\n";
	cout<<"Select the Type of Instance: \n\t\t\t0 for ordinary instance \n\t\t\t1 for Hard instance \n\t\t\t-1 for exit:		";
	cin>>d_selection;
	cout<<"Select the algorithm you want to use: \n\t\t\t0 for BBd \n\t\t\t1 for BBb \n\t\t\tOther for BBb2 :		";
	cin>>a_selection;
	if(a_selection == 0)
	{
		BBd_ = true;
		strcpy(Algorithm,"BBd");
	}
	else if(a_selection == 1)
	{
		BBb_ = true;
		strcpy(Algorithm,"BBb");
	}
	else
	{
		BBb2_= true;
		strcpy(Algorithm,"BBb2");
	}

	
		
	
	
	
	if(d_selection == -1)
		cout<<"Thanks for using this implementation\n";

	else if(d_selection == 1)
	{
		cout<<"For Hard instance, Please Enter the number of restriction sites n : ";
		cin>>n;
		problem_size = n;

		strcpy(fileName_output,"pdp_output_Zhang.txt");		
		writefile_time.open(fileName_time,ios::app);
		writefile_time<<"\n\n\t\tn\tm\tAlgorithm\tRunningTime\n";
		writefile_time.close();
		rand_solution = MakeHardInstance(n);
		rand_instance = ComputeDistance(rand_solution);
		print_vector(rand_solution);
		writefile_output.open(fileName_output,ios::app);
		writefile_output<<"\nn=	\t"<<n<<"\nOriginal Solution	=	";
		writefile_output.close();
		print_vector_in_file(rand_solution,fileName_output);
		writefile_output.open(fileName_output,ios::app);
		writefile_output<<"\nOriginal Problem	=	";
		writefile_output.close();
		print_vector_in_file(rand_instance,fileName_output);
		writefile_output.open(fileName_output,ios::app);
		writefile_output.close();
		//main function
		int start_s=clock();
		partialDigestProblem(rand_instance);
		int stop_s=clock();
		cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) <<"  Seconds"<< endl;
		runningTime = (stop_s-start_s)/double(CLOCKS_PER_SEC);
		writefile_time.open(fileName_time,ios::app);
		writefile_time<<"\t\t"<<rand_solution.size()<<"\t1000\t"<<Algorithm<<"\t"<<runningTime<<"\n";
		writefile_time.close();		

			
			
		
	}
	else if(d_selection == 0)
	{
		strcpy(fileName_output,"pdp_output_Random.txt");
		cout<<"For Ordinary instance Please Enter the number of restriction sites n : ";
		cin>>n;
		cout<<"Please Enter the average distance between resriction sites m : ";
		cin>>m;
		problem_size = n;

		writefile_time.open(fileName_time,ios::app);
		writefile_time<<"\n\n\t\tn\tm\tAlgorithm\tRunningTime\n";
		writefile_time.close();
		rand_solution = MakeRandomInstance(n, m);
		rand_instance = ComputeDistance(rand_solution);
		cout<<"\nSolution = "<<"	=	";
		print_vector(rand_solution);
		writefile_output.open(fileName_output,ios::app);
		writefile_output<<"\nn=	\t"<<n<<"m=	\t"<<m<<"\nOriginal Solution	=	";
		writefile_output.close();
		print_vector_in_file(rand_solution,fileName_output);
		writefile_output.open(fileName_output,ios::app);
		writefile_output<<"\nOriginal Problem	=	";
		writefile_output.close();
		print_vector_in_file(rand_instance,fileName_output);
		writefile_output.close();
		solutions.clear();
		int start_s=clock();
		//main function
		partialDigestProblem(rand_instance);
		int stop_s=clock();
		cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) <<"  Seconds"<< endl;
		runningTime = (stop_s-start_s)/double(CLOCKS_PER_SEC);
		writefile_time.open(fileName_time,ios::app);
		writefile_time<<"\n\n\t\t"<<rand_solution.size()<<"\t"<<m<<"\t"<<Algorithm<<"\t"<<runningTime<<"\n";
		writefile_time.close();
				
			
		

		
				
	}
	
	else
		cout<<"In correct choice Please run the program again and select correct choice\n";
}	
	else
	{
		
		BBb2_ = true;
		strcpy(Algorithm,"BBb2");
		int N=0;
		strcpy(fileName_output,"pdp_output_specific.txt");	
		vector< vector<int> > temp;
		ReadSet(argv[1], rand_instance, N);
		int n;
		n = int((1+sqrt((double)1+8*N))/2.0);
		cout<<"\n n = "<<n<<"\n";
		problem_size = n;
		writefile_time.open(fileName_time,ios::app);
		writefile_time<<"\n\n\t\tn\tAlgorithm\tRunningTime\n";
		writefile_time.close();
		sort(rand_instance.begin(),rand_instance.end());
		writefile_output.open(fileName_output,ios::app);
		writefile_output<<"\nOriginal Problem	=	";
		writefile_output.close();
		print_vector_in_file(rand_instance,fileName_output);
		writefile_output.open(fileName_output,ios::app);
		writefile_output.close();
		int start_s=clock();
		//main function
		partialDigestProblem(rand_instance);
		int stop_s=clock();
		cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) <<"  Seconds"<< endl;
		runningTime = (stop_s-start_s)/double(CLOCKS_PER_SEC);

		//print_solutions_in_file (temp,fileName_output);
		cout<<"Time  "<<runningTime<<"\n";
		writefile_time.open(fileName_time,ios::app);
		writefile_time<<"\n\n\t\t"<<n<<"\tBBb2\t"<<Algorithm<<"\t"<<runningTime<<"\n";
		writefile_time.close();

	}

	int xx;
	
	cin>>xx;
	return 0;
}