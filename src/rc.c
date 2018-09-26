void rc(char *str){
	char tmp;
	long n,i;

	n = strlen(str);
	i=0;
	while(i<=ceil(n/2)){
		tmp=str[i];
		switch(str[n-i-1]){
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
			str[i]=str[n-i-1];
		}
		switch(tmp){
		case 'A':
			str[n-i-1]='T';
			break;
		case 'C':
			str[n-i-1]='G';
			break;
		case 'G':
			str[n-i-1]='C';
			break;
		case 'T':
			str[n-i-1]='A';
			break;
		default:
			str[n-i-1]=tmp;
		}
		i++;
	}
}
