#include <math.h>
#include <string.h>
#include <stdio.h>

void rc(char *str){
	char tmp;
	unsigned long n,i,j;
	unsigned long middle;

	n = strlen(str);
	middle=floor(n/2); //If n is divisible by 2  this is leftmost index of the right side eg: 01|23 4/2=2
	                  // if n is not divisble by 01|2|34 5/2 =2.5 rounds to 2
	                  // In all cases we don't need to consider it.
	i=0;              // Left end
	j=n-1;             // right end
	//printf("String length is %ld and middle is %ld\n",n,middle);
	while(i<middle){
//		printf("Swapping %c with %c at %ld and %ld\n",str[i],str[j],i,j);
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
