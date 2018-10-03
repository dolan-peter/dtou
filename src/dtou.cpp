#include <Rcpp.h>
using namespace Rcpp;

extern void rc(char *str);

char *S1;        // Global character array containing String imported from R
long *I1;        // Global Stack
long *dtou;      // Global array that will hold the results
long i,j,k,next; //We make them global so space isn't wasted inside the recursion routine.
long d;          // A depth cutoff [only used by depth-limited version]
long l, m;       // globals used by the stack-size variant
long _stringLength=0; // A global used to keep track of the size (change based on whether RC is being analyzed or not)
bool verbose;

//Global Initialization:
List convertRStringsandGlobalSetup(std::vector<std::string> RString, bool parseRC, int depth){
	/*
	* Converts list of Rstrings into single X-delimted string stored in global S1
	* Also prepares an R List to hold the results of the calculation
	*      the int depth is only used by some of the function calls
	* SideEffects: Intializes String S1
	*            : Sets n to the appropriate value
	*            :
	*/


	long unsigned n,i,j,k,tmp,n1;
	long unsigned n_strings=RString.size();

	List results=List::create();

	n=0;
	for(j=0;j<n_strings;j++){
		tmp = RString[j].size()+1;
		n+=tmp; // We need the +1 to make space for the X's that will delimit the strings
		results.push_back(NumericVector(tmp-1)); // Add an empty Numeric Vector for the results
	}
	n1=parseRC?2*n:n; // Prepare allocation sizes for rc vs plain vanilla

	Rcout << "Preparing to calculate dtou for " << n << " nucleotides and boundaries...\n";

	S1=new char[n1+1]; //S1 is a global
	if(S1==0){stop("Unable to allocate space for character string\n");}

	k=0; // It's already 0... but it's nice to be explicit
	for(j=0;j< n_strings;j++){
		for(i=0;i<RString[j].size();i++){
			S1[k]=RString[j][i];
			k++;
		}
		S1[k++]='X'; // Add the X at the end
	}
	Rcout << "\tProcessed " << k << " characters/boundaries.  Expected " << n <<".\n";
	S1[n] = 0; // If we're not parsing RC then we're ready to go... otherwise this will be overwritten:
  //Rcout << "Before duplication " <<S1<<"\n";
	if(parseRC){
		for(i=0;i<n;i++){ // duplicate string
			S1[n+i]=S1[i];
		}
		S1[2*n]=0;        // Add a safety margin.
//		Rcout << "After duplication " <<S1<<"\n";
		rc(S1+n-1);       // This way we include the final 'X'
	}
//	Rcout << "After rc " <<S1<<"\n";
	I1=new long[n1+1];
	dtou = new long[n1+1]; // We should not need to initialize these-- that should be handled automatically during processing

	// Initialize Array encoding virtual tree of stacks:
	for(i=0;i<n1;i++){I1[i]=i+1;}
	I1[n1]=-1;

	d=depth;
	_stringLength=n1;
	return(results);
}

void rc(char *str){
	char tmp;
	unsigned long n,i,j;

	n = strlen(str);

	i=0;              // Left end
	j=n-1;             // right end
	while(i<=j){
		tmp=str[i];     // store the left end
		switch(str[j]){ // Now look at the right
		case 'A':
			str[i]='T';
			break;
		case 'C':
			str[i]='G';
			break;
		case 'G':
			str[i]='C';
			break;
		case 'T':
			str[i]='A';
			break;
		default:
			str[i]=str[j]; // Stays the same if it's not ACGT
		}
		switch(tmp){
		case 'A':
			str[j]='T';
			break;
		case 'C':
			str[j]='G';
			break;
		case 'G':
			str[j]='C';
			break;
		case 'T':
			str[j]='A';
			break;
		default:
			str[j]=tmp;
		}
		i++;j--;
	}
}

// The workhorses-- there is a LOT of code redundancy-- but in the interests of efficiency I am not adding conditionals
void recurse(long Ip1, long depth){
	long A1=-1, C1=-1, G1=-1, T1=-1;
	if(Ip1==-1){
		return;}
	if(I1[Ip1]==-1){ // Only one entry exists in the stack
		dtou[Ip1]=depth;
		return;
	}

	i=Ip1;
	while(i!=-1){
		next=I1[i];
		switch(S1[i+depth]){
		case 'A':
			I1[i]=A1;
			A1=i;
			break;
		case 'C':
			I1[i]=C1;
			C1=i;
			break;
		case 'G':
			I1[i]=G1;
			G1=i;
			break;
		case 'T':
			I1[i]=T1;
			T1=i;
			break;
		case 'N':
		case 'X':
			I1[i]=-1;
			dtou[i]=depth;
			break;
		}
		i=next;
	}
	recurse(A1,depth+1);
	recurse(C1,depth+1);
	recurse(G1,depth+1);
	recurse(T1,depth+1);
	return;
}

void depthLimitedRecurse(long Ip1, long depth){
	long A1=-1, C1=-1, G1=-1, T1=-1;

	if(depth>d){return;}
	if(Ip1==-1){return;}

	if(I1[Ip1]==-1){ // Only one entry exists in the stack
		dtou[Ip1]=depth;
		return;
	}

	i=Ip1;
	while(i!=-1){
		dtou[i]=depth; //Assign all the values in the current stack to the current depth
		next=I1[i];
		switch(S1[i+depth]){
		case 'A':
			I1[i]=A1;
			A1=i;
			break;
		case 'C':
			I1[i]=C1;
			C1=i;
			break;
		case 'G':
			I1[i]=G1;
			G1=i;
			break;
		case 'T':
			I1[i]=T1;
			T1=i;
			break;
		case 'N':
		case 'X':
			I1[i]=-1;
			dtou[i]=depth;
			break;
		}
		i=next;
	}

	depthLimitedRecurse(A1,depth+1);
	depthLimitedRecurse(C1,depth+1);
	depthLimitedRecurse(G1,depth+1);
	depthLimitedRecurse(T1,depth+1);
	return;
}

void recurseWithSize(long Ip1, long depth,long stackSize){
	long A1=-1, C1=-1, G1=-1, T1=-1;
	long ASize=0,CSize=0,GSize=0,TSize=0;

	if(Ip1==-1){
		return;}
	if(I1[Ip1]==-1){ // Only one entry exists in the stack
		dtou[Ip1]=depth;
		return;
	}
	if(stackSize==2){
		l = Ip1+depth;
		m = I1[Ip1]+depth;

		while(S1[l]==S1[m]){
			l++;m++;
		} ; // Run until they are different... assumption is that we will NOT cross 'X's
		// That might not be warranted-- I'll try it without the safeguard first
		// And then see what kind of performance hit I get if I scan for 'X'
		dtou[Ip1]=l-Ip1+1;
		if(S1[l]=='X'){dtou[Ip1]--;}

		dtou[I1[Ip1]]=l-Ip1+1;
		if(S1[m]=='X'){dtou[I1[Ip1]]--;}

		I1[Ip1]=-1;               // Close down this stack... other value should already be -1
		return;
	}

	i=Ip1;
	while(i!=-1){
		next=I1[i];
		switch(S1[i+depth]){
		case 'A':
			I1[i]=A1;
			A1=i;
			ASize++;
			break;
		case 'C':
			I1[i]=C1;
			C1=i;
			CSize++;
			break;
		case 'G':
			I1[i]=G1;
			G1=i;
			GSize++;
			break;
		case 'T':
			I1[i]=T1;
			T1=i;
			TSize++;
			break;
		case 'N':
		case 'X':
			I1[i]=-1;
			dtou[i]=depth;
			break;
		}
		i=next;
	}
	recurseWithSize(A1,depth+1,ASize);
	recurseWithSize(C1,depth+1,CSize);
	recurseWithSize(G1,depth+1,GSize);
	recurseWithSize(T1,depth+1,TSize);
	return;
}

void recurseWithSizeB(long Ip1, long depth,long stackSize){
	long A1=-1, C1=-1, G1=-1, T1=-1;
	long ASize=0,CSize=0,GSize=0,TSize=0;

	if(Ip1==-1){
		return;}
	if(I1[Ip1]==-1){ // Only one entry exists in the stack
		dtou[Ip1]=(dtou[Ip1]>depth)?dtou[Ip1]:depth;
		return;
	}
	if(stackSize==2){
		l = Ip1+depth;
		m = I1[Ip1]+depth;

		if(dtou[Ip1]>depth && dtou[I1[Ip1]]>depth){ // If we've already parsed these regions and they're longer than our current depth than this run is superfluous
			return;
		}
		while(S1[l]==S1[m]){
			l++;m++;
		} ; // Run until they are different... assumption is that we will NOT cross 'X's
		// That might not be warranted-- I'll try it without the safeguard first
		// And then see what kind of performance hit I get if I scan for 'X'

		// Now we need to set values for dtou
		//	dtou[I1[Ip1]]=l-Ip1+1;
		// dtou[Ip1]=l-Ip1+1;
		m=0;
		while(l>Ip1+depth+1){ // Don't go too deep with the over-writing
			dtou[I1[Ip1]+m]=(dtou[I1[Ip1]+m]>l-Ip1+1)?(dtou[I1[Ip1]+m]):(l-Ip1+1);
			dtou[Ip1+m]=(dtou[Ip1+m]>l-Ip1+1)?(dtou[Ip1+m]):(l-Ip1+1);
			l--;
			m++;
		}
		if(S1[l]=='X'){dtou[Ip1]--;}

		dtou[I1[Ip1]]=l-Ip1+1;
		if(S1[m]=='X'){dtou[I1[Ip1]]--;}

		I1[Ip1]=-1;               // Close down this stack... other value should already be -1
		return;
	}

	i=Ip1;
	while(i!=-1){
		next=I1[i];
		switch(S1[i+depth]){
		case 'A':
			I1[i]=A1;
			A1=i;
			ASize++;
			break;
		case 'C':
			I1[i]=C1;
			C1=i;
			CSize++;
			break;
		case 'G':
			I1[i]=G1;
			G1=i;
			GSize++;
			break;
		case 'T':
			I1[i]=T1;
			T1=i;
			TSize++;
			break;
		case 'N':
		case 'X':
			I1[i]=-1;
			dtou[i]=depth;
			break;
		}
		i=next;
	}
	recurseWithSizeB(A1,depth+1,ASize);
	recurseWithSizeB(C1,depth+1,CSize);
	recurseWithSizeB(G1,depth+1,GSize);
	recurseWithSizeB(T1,depth+1,TSize);
	return;
}

void depthLimitedRecurseWithSize(long Ip1, long depth,long stackSize){
	long A1=-1, C1=-1, G1=-1, T1=-1;
	long ASize=0,CSize=0,GSize=0,TSize=0;

	if(depth>d){return;}
	if(Ip1==-1){return;}
	if(I1[Ip1]==-1){ // Only one entry exists in the stack
		dtou[Ip1]=depth;
		return;
	}
	if(stackSize==2){
		l = Ip1+depth;
		m = I1[Ip1]+depth;

		while(l-Ip1 <= d && S1[l]==S1[m]){
			l++;m++;
		} ; // Run until they are different... assumption is that we will NOT cross 'X's
		// That might not be warranted-- I'll try it without the safeguard first
		// And then see what kind of performance hit I get if I scan for 'X'
		dtou[Ip1]=l-Ip1+1;
		if(S1[l]=='X'){dtou[Ip1]--;}

		dtou[I1[Ip1]]=l-Ip1+1;
		if(S1[m]=='X'){dtou[I1[Ip1]]--;}

		I1[Ip1]=-1;               // Close down this stack... other value should already be -1
		return;
	}

	i=Ip1;
	while(i!=-1){
		dtou[i]=depth; //Assign all the values in the current stack to the current depth... not the most efficient solution
		next=I1[i];
		switch(S1[i+depth]){
		case 'A':
			I1[i]=A1;
			A1=i;
			ASize++;
			break;
		case 'C':
			I1[i]=C1;
			C1=i;
			CSize++;
			break;
		case 'G':
			I1[i]=G1;
			G1=i;
			GSize++;
			break;
		case 'T':
			I1[i]=T1;
			T1=i;
			TSize++;
			break;
		case 'N':
		case 'X':
			I1[i]=-1;
			dtou[i]=depth;
			break;
		}
		i=next;
	}
	depthLimitedRecurseWithSize(A1,depth+1,ASize);
	depthLimitedRecurseWithSize(C1,depth+1,CSize);
	depthLimitedRecurseWithSize(G1,depth+1,GSize);
	depthLimitedRecurseWithSize(T1,depth+1,TSize);
	return;
}

void updateResults(std::vector<std::string> RString,List results){
	long unsigned i,j,tmp;
	long unsigned n=0;
	long unsigned n_strings=RString.size();

	for(j=0;j<n_strings;j++){
		tmp=RString[j].size();
		for(i=0; i<tmp;i++){
			((NumericVector) results[j])[i]=dtou[n];
			n++;
		}
		n++; // We do NOT record the results of the 'X'... but we do need to pass-over the result

	}
}


// [[Rcpp::export]]
List c_dtou(std::vector<std::string> RString,bool rc){ //TODO:  Ensure that RString is using ASCII encoding
	  List results=convertRStringsandGlobalSetup(RString,rc,0); // This version doesn't use-depth-limiting

		Rcout << "\tBeginning scan...";
		recurse(0,0);
		Rcout <<"done\n";

		updateResults(RString,results); //Prepares results for being returned-- needs RString to get the sizes correct

	  return results;
}

// [[Rcpp::export]]
List c_dtouDepthLimit(std::vector<std::string> RString,bool rc,long depth){ //TODO:  Ensure that RString is using ASCII encoding
	List results=convertRStringsandGlobalSetup(RString,rc,depth); // This version use-depth-limiting

	Rcout << "\tBeginning scan...";
	depthLimitedRecurse(0,0);
	Rcout <<"done\n";

	updateResults(RString,results); //Prepares results for being returned-- needs RString to get the sizes correct

	return results;
}

// [[Rcpp::export]]
List c_dtouS2(std::vector<std::string> RString,bool rc){ //TODO:  Ensure that RString is using ASCII encoding
	List results=convertRStringsandGlobalSetup(RString,rc,0); // This version doesn't use-depth-limiting

	Rcout << "\tBeginning scan...";
	recurseWithSize(0,0,_stringLength+1);
	Rcout <<"done\n";

	updateResults(RString,results); //Prepares results for being returned-- needs RString to get the sizes correct

	return results;
}

// [[Rcpp::export]]
List c_dtouS2B(std::vector<std::string> RString,bool rc){ //TODO:  Ensure that RString is using ASCII encoding
	List results=convertRStringsandGlobalSetup(RString,rc,0); // This version doesn't use-depth-limiting

	Rcout << "\tBeginning scan...";
	recurseWithSizeB(0,0,_stringLength+1);
	Rcout <<"done\n";

	updateResults(RString,results); //Prepares results for being returned-- needs RString to get the sizes correct

	return results;
}


// Change NumericVector to list
// [[Rcpp::export]]
List c_dtouS2DepthLimit(std::vector<std::string> RString, bool rc,long depth){ //TODO:  Ensure that RString is using ASCII encoding
	List results=convertRStringsandGlobalSetup(RString,rc,depth); // This version is depth-limiting

	Rcout << "\tBeginning scan...";
	depthLimitedRecurseWithSize(0,0,_stringLength+1);
	Rcout <<"done\n";

	updateResults(RString,results); //Prepares results for being returned-- needs RString to get the sizes correct

	return results;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
c_dtou(c("AAAAACCCGACTGGGCTCA","ACCT"),TRUE)
c_dtou("AAAAACCCGACTGGGCTCAXACCT",TRUE)
c_dtou("AAAAACCCGACTGGGCTCAXACCT",FALSE)
c_dtouDepthLimit("AAAAACCCGACTGGGCTCAX",TRUE,10)
c_dtouDepthLimit(c("AAAAACCCGACTGGGCTCA","ACCT"),TRUE,2)
c_dtouDepthLimit("AAAAACCCGACTGGGCTCAXACCT",FALSE,4)
c_dtouDepthLimit("AAAAACCCGACTGGGCTCAXACCT",TRUE,4)
c_dtouS2("AAAAACCCGACTGGGCTCA",TRUE)
c_dtouS2("AAAAACCCGACTGGGCTCA",FALSE)
c_dtouS2(c("AAAAACCCGACTGGGCTCA","ACCT"),TRUE)
c_dtouS2DepthLimit("AAAAACCCGACTGGGCTCAX",TRUE,10)
c_dtouS2DepthLimit(c("AAAAACCCGACTGGGCTCA","ACCT"),TRUE,2)
c_dtouS2DepthLimit("AAAAACCCGACTGGGCTCAXACCT",FALSE,4)
*/
