#include <Rcpp.h>
using namespace Rcpp;

extern void rc(char *str);

char *str;        // Global character array containing String imported from R
long *Indices;        // Global Stack
long counts[101];  // Global array that will hold the results
long _i,_k,_next; //We make them global so space isn't wasted inside the recursion routine.
long _d;          // A depth cutoff [only used by depth-limited version]
long __stringLength=0; // A global used to keep track of the size (change based on whether RC is being analyzed or not)

//Global Initialization:
void VconvertRStringsandGlobalSetup(std::vector<std::string> RString, bool parseRC, int depth){
	long unsigned n,i,j,k,tmp,n1;
	long unsigned n_strings=RString.size();

	n=0;
	for(j=0;j<n_strings;j++){
		tmp = RString[j].size()+1;
		n+=tmp; // We need the +1 to make space for the X's that will delimit the strings
	}
	n1=parseRC?2*n:n; // Prepare allocation sizes for rc vs plain vanilla

	Rcout << "Preparing to calculate dtou for " << n << " nucleotides and boundaries...\n";

	str=new char[n1+1]; //str is a global
	if(str==0){stop("Unable to allocate space for character string\n");}

	k=0; // It's already 0... but it's nice to be explicit
	for(j=0;j< n_strings;j++){
		for(i=0;i<RString[j].size();i++){
			str[k]=RString[j][i];
			k++;
		}
		str[k++]='X'; // Add the X at the end
	}
	str[n] = 0; // If we're not parsing RC then we're ready to go... otherwise this will be overwritten:
	if(parseRC){
		for(i=0;i<n;i++){ // duplicate string
			str[n+i]=str[i];
		}
		str[2*n]=0;        // Add a safety margin.
		rc(str+n-1);       // This way we include the final 'X'
	}
	Indices=new long[n1+1];

	// Initialize Array encoding virtual tree of stacks:
	for(i=0;i<n1;i++){Indices[i]=i+1;}
	Indices[n1]=-1;

	_d=depth;
	__stringLength=n1;

	for(i=0;i<=100;i++){counts[i]=0;}
}

void depthLimitedVignetteRecurse(long Ip1, long depth,long stackSize){
	long A1=-1, C1=-1, G1=-1, T1=-1;
	long ASize=0,CSize=0,GSize=0,TSize=0;

	if(depth>_d){return;}
	if(Ip1==-1){return;}
	counts[depth]+=stackSize;
	if(Indices[Ip1]==-1){ // Only one entry exists in the stack
		return;
	}
	_i=Ip1;
	while(_i!=-1){
		_next=Indices[_i];
		switch(str[_i+depth]){
		case 'A':
			Indices[_i]=A1;
			A1=_i;
			ASize++;
			break;
		case 'C':
			Indices[_i]=C1;
			C1=_i;
			CSize++;
			break;
		case 'G':
			Indices[_i]=G1;
			G1=_i;
			GSize++;
			break;
		case 'T':
			Indices[_i]=T1;
			T1=_i;
			TSize++;
			break;
		case 'N':
		case 'X':
			Indices[_i]=-1;
			break;
		}
		_i=_next;
	}
	depthLimitedVignetteRecurse(A1,depth+1,ASize);
	depthLimitedVignetteRecurse(C1,depth+1,CSize);
	depthLimitedVignetteRecurse(G1,depth+1,GSize);
	depthLimitedVignetteRecurse(T1,depth+1,TSize);
	return;
}

// [[Rcpp::export]]
NumericVector c_vignetteExample(std::vector<std::string> RString,bool rc){
	long int i;
	VconvertRStringsandGlobalSetup(RString,rc,100); // This version is depth-limiting up to depth 100

	Rcout << "\tBeginning Recursion...";
	depthLimitedVignetteRecurse(0,0,__stringLength+1);
	Rcout << "done\n";
  Rcout << "\tPreparing output...";
	NumericVector outputR(101);
	for(i=0; i<=100;i++){
		outputR[i]=counts[i];
	}
	Rcout << "done\n";
	return outputR;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
c_vignetteExample("AACGCCCGCGCTCCCGCCGCCCG",FALSE)
*/
