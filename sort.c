#include <stdio.h>
#include <stdlib.h>
#define N 196
#define numOfFeature 196

int main(){
	double myData[196],data[196][196];
	int i,j,z,l,next;
	int dist;
	int group=1;
	int member;
	
	FILE *fpw,*fpr;
	char fnw[40],fnr[40];
	int saveOrder[196];
    double tmp[196][196];
    int nextOrder;
  

	for(l=0;l<46;l++){
		sprintf(fnr,"./vectorData/value%02d.txt",l+1);
		sprintf(fnw,"./sortedData/sortedValue%02d.txt",l+1);
		printf("open:%s out:%s\n",fnr,fnw);
		fpr=fopen(fnr,"r");
		fpw=fopen(fnw,"w");
		if(fpr==NULL){
			printf("cannot open a file[read]\n");
			exit(1);
		}
		if(fpw==NULL){
			printf("cannot open a file[write]\n");
			exit(1);
		}
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				fscanf(fpr,"%lf",&data[i][j]);
			}
		}
		for(i=0;i<N;i++){
			myData[i]=data[i][i];
			saveOrder[i]=i;
		}
		while(group*3+1<N)group=group*3+1;
		while(group>=1){
			member=N/group;
			for(z=0;z<group;z++){
				for(i=1;i<=member-1;i++){
					dist=z+group*i;
					next=myData[dist];
	                nextOrder=saveOrder[dist];
					for(j=dist;j>=z+group && myData[j-group] < next;j=j-group){
						myData[j] = myData[j-group];
                    	saveOrder[j] = saveOrder[j-group];
					}
					myData[j]=next;
	                saveOrder[j]=nextOrder;
				}
			}
			group=group/3;
		}
		for(i=0;i<N;i++){
			fprintf(fpw,"%lf\n",myData[i]);
		}
		fclose(fpw);
		fclose(fpr);


		sprintf(fnr,"./vectorData/vector%02d.txt",l+1);
		sprintf(fnw,"./sortedData/sortedVector%02d.txt",l+1);
		printf("open:%s out:%s\n",fnr,fnw);
		fpr=fopen(fnr,"r");
		fpw=fopen(fnw,"w");
		if(fpr==NULL){
			printf("cannot open a file[read]\n");
			exit(1);
		}
		if(fpw==NULL){
			printf("cannot open a file[write]\n");
			exit(1);
		}
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				fscanf(fpr,"%lf",&data[i][j]);
			}
		}

		for(i=0;i<numOfFeature;i++){
	        for(j=0;j<numOfFeature;j++){
    	        tmp[i][j]=data[saveOrder[i]][j];
    	    }
	    }
		for(i=0;i<numOfFeature;i++){
			for(j=0;j<numOfFeature;j++){
				data[i][j]=tmp[i][j];
			}
		}
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				fprintf(fpw,"%lf ",data[i][j]);
			}
			fprintf(fpw,"\n");
		}
		fclose(fpw);
		fclose(fpr);
	}
	return 0;
}



