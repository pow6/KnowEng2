//4J02 s15015 �r�����i
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
    char fnReadFormat[12]="c",fnWriteFormat[12]="mean";
    char fnRead[24],fnWrite[24];
    int numOfData=180,i,j;
    double (*data)[numOfFeature],average[numOfFeature];
    data = malloc(sizeof(double)*numOfData*numOfFeature);
    for(i=0;i<46;i++){
        sprintf(fnRead,"%s%02d.txt",fnReadFormat,i+1);
        printf("Read[%s]\n",fnRead);
        readData(numOfData,data,fnRead);
        average_calcFeature(data,average,numOfData);
        sprintf(fnWrite,"%s%02d.txt",fnWriteFormat,i+1);
        printf("Write[%s]\n",fnWrite);
        writeData(numOfData,average,fnWrite);
    }
}

//�t�@�C������f�[�^��ǂݎ��
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

//�t�@�C���ɔz�����������
void writeData(int numOfData,double average[],char fileName[]){
    FILE *write;
    int i;
    fileWrite(fileName,write);
    for(i=0;i<numOfData;i++){
        fprintf(write,"%lf\n",average[i]);
    }
    fclose(write);
}

//�P�̗v�f���̕��ϓ����ʂ��v�Z
double oneElement_calcFeature(double data[][numOfFeature],int numOfData,int target){
    int i;
    double result=0;
    for(i=0;i<numOfData;i++){
        result+=data[i][target];
    }
    return result/numOfData;
}

//�w�萔�̕������̕��ϓ����ʂ��v�Z
void average_calcFeature(double data[][numOfFeature],double average[],int numOfData){
    int i;
    for(i=0;i<numOfFeature;i++){
        average[i]=oneElement_calcFeature(data,numOfData,i);
    }
}

