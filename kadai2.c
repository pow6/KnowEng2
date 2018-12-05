//4J02 s15015 �r�����i
#include <stdio.h>
#include <stdlib.h>

#define numOfFeature 196
#define numOfData 180
#define fileRead(fileName,fileStream) fileStream=fopen(fileName,"r");if(fileStream==NULL){printf("cannat read file[%s]",fileName);exit(1);}
#define fileWrite(fileName,fileStream) fileStream=fopen(fileName,"w");if(fileStream==NULL){printf("cannat write file[%s]",fileName);exit(1);}

void readData(double data[][numOfFeature],char fileName[]);
void writeData(double average[],char fileName[]);
void writeDataTwoDim(double coriance[][numOfFeature],char fileName[]);
double oneElement_calcFeature(double data[][numOfFeature],int target);
void average_calcFeature(double data[][numOfFeature],double average[]);
void calcCovariance(double data[][numOfFeature],double average[],double covariance[][numOfFeature]);

void main(){
    char fnReadFormat[12]="c",fnWriteFormat[12]="sigma";
    char fnRead[24],fnWrite[24];
    int i,j;
    double (*data)[numOfFeature],average[numOfFeature],(*covariance)[numOfFeature];
    data = malloc(sizeof(double)*numOfData*numOfFeature);
    covariance = malloc(sizeof(double)*numOfFeature*numOfFeature);
    for(i=0;i<46;i++){
        sprintf(fnRead,"./originData/%s%02d.txt",fnReadFormat,i+1);
        readData(data,fnRead);
        average_calcFeature(data,average);
        calcCovariance(data,average,covariance);
        sprintf(fnWrite,"./sigmaData/%s%02d.txt",fnWriteFormat,i+1);
        writeDataTwoDim(covariance,fnWrite);
    }
}

//�t�@�C������f�[�^��ǂݎ��
void readData(double data[][numOfFeature],char fileName[]){
    FILE *read;
    int i,j;
    fileRead(fileName,read);
    printf("Read[%s]\n",fileName);
    for(j=0;j<numOfData;j++){
        for(i=0;i<numOfFeature;i++){
            fscanf(read,"%lf",&data[j][i]);
        }
    }
    fclose(read);
}

//�t�@�C���ɔz�����������
void writeData(double average[],char fileName[]){
    FILE *write;
    int i;
    fileWrite(fileName,write);
    printf("Write[%s]\n",fileName);
    for(i=0;i<numOfData;i++){
        fprintf(write,"%lf\n",average[i]);
    }
    fclose(write);
}

//�t�@�C���ɂQ�����z�����������
void writeDataTwoDim(double coriance[][numOfFeature],char fileName[])
{
    FILE *write;
    int i,j;
    fileWrite(fileName,write);
    printf("Write[%s]\n",fileName);
    for(i=0;i<numOfFeature;i++){
        for(j=0;j<numOfFeature;j++){
            fprintf(write,"%lf ",coriance[i][j]);
        }
        fprintf(write,"\n");
    }
    fclose(write);
}

//�P�̗v�f���̕��ϓ����ʂ��v�Z
double oneElement_calcFeature(double data[][numOfFeature],int target){
    int i;
    double result=0;
    for(i=0;i<numOfData;i++){
        result+=data[i][target];
    }
    return result/numOfData;
}

//�w�萔�̕������̕��ϓ����ʂ��v�Z(���)
void average_calcFeature(double data[][numOfFeature],double average[]){
    int i;
    for(i=0;i<numOfFeature;i++){
        average[i]=oneElement_calcFeature(data,i);
    }
}

//�����U���v�Z����i�e�v�f���ƂɌv�Z�j
void calcCovariance(double data[][numOfFeature],double average[],double covariance[][numOfFeature])
{
    int k,j,i;
    for(i=0;i<numOfFeature;i++){
        for(j=0;j<numOfFeature;j++){
            covariance[i][j]=0;
        }
    }
    for(i=0;i<numOfFeature;i++){
        for(j=0;j<numOfFeature;j++){
            for(k=0;k<numOfData;k++){
                covariance[i][j]+=data[k][i]*data[k][j];
            }
            covariance[i][j]=covariance[i][j]/numOfData-average[i]*average[j];
        }
    }
}

