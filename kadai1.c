//4J02 s15015 池口恭司
#include <stdio.h>
#include <stdlib.h>

#define numOfFeature 196
#define fileRead(fileName,fileStream) fileStream=fopen(fileName,"r");if(fileStream==NULL){printf("cannat read file[%s]",fileName);exit(1);}
#define fileWrite(fileName,fileStream) fileStream=fopen(fileName,"w");if(fileStream==NULL){printf("cannat write file[%s]",fileName);exit(1);}

void readData(int numOfData,double data[][numOfFeature],char fileName[]);
void writeData(int numOfData,double average[],char fileName[]);
double oneElement_calcFeature(double data[][numOfFeature],int numOfData,int target);
void average_calcFeature(double data[][numOfFeature],double average[],int numOfData);

void main(){
    char fnReadFormat[12],fnWriteFormat[12];
    char fnRead[24],fnWrite[24];
    int numOfData,i,j;
    double (*data)[numOfFeature],average[numOfFeature];
    printf("読み込むデータ数：\n");
    scanf("%d",&numOfData);
    printf("入力ファイルフォーマット（数字前の先頭の文字を入力/c01.txt~c46.txtなら、cを入力）：");
    scanf("%s",fnReadFormat);
    printf("出力ファイルフォーマット（数字前の先頭の文字を入力/mean01.txt~mean46.txtなら、meanを入力）：");
    scanf("%s",fnWriteFormat);
    data = malloc(sizeof(double)*numOfData*numOfFeature);
    for(i=0;i<46;i++){
        //sprintf(fnRead,"%s%02d.txt",fnReadFormat,i+1);
        sprintf(fnRead,"./%s%02d.txt",fnReadFormat,i+1);
        printf("Read[%s]\n",fnRead);
        readData(numOfData,data,fnRead);
        average_calcFeature(data,average,numOfData);
        //sprintf(fnWrite,"%s%02d.txt",fnWriteFormat,i+1);
        sprintf(fnWrite,"./meanData/%s%02d.txt",fnWriteFormat,i+1);
        printf("Write[%s]\n",fnWrite);
        writeData(numOfData,average,fnWrite);
    }
    
}

//ファイルからデータを読み取る
void readData(int numOfData,double data[][numOfFeature],char fileName[]){
    FILE *read;
    int i,j;
    fileRead(fileName,read);
    for(j=0;j<numOfData;j++){
        for(i=0;i<numOfFeature;i++){
            fscanf(read,"%lf",&data[j][i]);
        }
    }
    fclose(read);
}

//ファイルに配列を書き込む
void writeData(int numOfData,double average[],char fileName[]){
    FILE *write;
    int i;
    fileWrite(fileName,write);
    for(i=0;i<numOfData;i++){
        fprintf(write,"%lf\n",average[i]);
    }
    fclose(write);
}

//１つの要素分の平均特徴量を計算
double oneElement_calcFeature(double data[][numOfFeature],int numOfData,int target){
    int i;
    double result=0;
    for(i=0;i<numOfData;i++){
        result+=data[i][target];
    }
    return result/numOfData;
}

//指定数の文字分の平均特徴量を計算
void average_calcFeature(double data[][numOfFeature],double average[],int numOfData){
    int i;
    for(i=0;i<numOfFeature;i++){
        average[i]=oneElement_calcFeature(data,numOfData,i);
    }
}

