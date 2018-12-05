//4J02 s15015 池口恭司
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define zero 0.01
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

void main()
{
    char fnReadFormat[12]="c",fnWriteFormat[12]="sigma";
    char fnRead[24],fnWrite[24];
    int i,j;
    double (*data)[numOfFeature],average[numOfFeature],(*covariance)[numOfFeature],(*eigenvalue)[numOfFeature];
    data = malloc(sizeof(double)*numOfData*numOfFeature);
    covariance = malloc(sizeof(double)*numOfData*numOfFeature);
    eigenvalue = malloc(sizeof(double)*numOfData*numOfFeature);
    for(i=0;i<46;i++){
        sprintf(fnRead,"./originData/%s%02d.txt",fnReadFormat,i+1);
        readData(data,fnRead);
        average_calcFeature(data,average);
        calcCovariance(data,average,covariance);
        sprintf(fnWrite,"./sigmaData/%s%02d.txt",fnWriteFormat,i+1);
        writeDataTwoDim(covariance,fnWrite);
    }
}

//ファイルからデータを読み取る
void readData(double data[][numOfFeature],char fileName[])
{
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

//ファイルに配列を書き込む
void writeData(double average[],char fileName[])
{
    FILE *write;
    int i;
    fileWrite(fileName,write);
    printf("Write[%s]\n",fileName);
    for(i=0;i<numOfData;i++){
        fprintf(write,"%lf\n",average[i]);
    }
    fclose(write);
}

//ファイルに２次元配列を書き込む
void writeDataTwoDim(double coriance[][numOfFeature],char fileName[])
{
    FILE *write;
    int i,j;
    fileWrite(fileName,write);
    printf("Write[%s]\n",fileName);
    for(i=0;i<numOfData;i++){
        for(j=0;j<numOfFeature;j++){
            fprintf(write,"%lf ",coriance[i][j]);
        }
        fprintf(write,"\n");
    }
    fclose(write);
}

//１つの要素分の平均特徴量を計算
double oneElement_calcFeature(double data[][numOfFeature],int target)
{
    int i;
    double result=0;
    for(i=0;i<numOfData;i++){
        result+=data[i][target];
    }
    return result/numOfData;
}

//指定数の文字分の平均特徴量を計算(代入)
void average_calcFeature(double data[][numOfFeature],double average[])
{
    int i;
    for(i=0;i<numOfFeature;i++){
        average[i]=oneElement_calcFeature(data,i);
    }
}

//共分散を計算する（各要素ごとに計算）
void calcCovariance(double data[][numOfFeature],double average[],double covariance[][numOfFeature])
{
    int k,j,i;
    for(i=0;i<numOfData;i++){
        for(j=0;j<numOfFeature;j++){
            covariance[i][j]=0;
        }
    }
    for(i=0;i<numOfData;i++){
        for(j=0;j<numOfFeature;j++){
            for(k=0;k<numOfData;k++){
                covariance[i][j]+=data[k][i]*data[k][j];
            }
            covariance[i][j]=covariance[i][j]/numOfData-average[i]*average[j];
        }
    }
}

//固有値・固有ベクトルの計算
//(196次元+1)個の点→平面が１つ決まる
//しかし、今回は1〜180個のデータのみ利用する
//本来であれば、197個ないと固有値・固有ベクトルは求まらないが、計算誤差のため、求まる
//誤差がある値なので、正規直交系はe1~e180までを利用することにする（基本的には）
//ただ、求まらない場合もあるため、求まらなくなった時点で計算を打ち切る必要がある
void calcEigenvalue(double covariance[][numOfFeature],double eigenvalue[][numOfFeature])
{
    int k,j,i;
    for(i=0;i<numOfData;i++){
        for(j=0;j<numOfFeature;j++){
            eigenvalue[i][j]=covariance[i][j];
        }
    }
    
}

//一つの文字あたりの固有値の計算処理について、計算を打ち切るなどの判定を行う
void oneWord_calcEigenvalue(double eigenvalue[][numOfFeature])
{
    
}

int oneWord_calcEigenvalueExe(double eigenvalue[][numOfFeature],int small,int big)  //(i,j)を(small,big)で表す
{
    int i,j;
    static int counter=0;
    double (*tmp)[numOfFeature];
    double valueOfSin,valueOfCos,theta,mulSinSS,mulSinBB,mulSinSB,mulCosSS,mulCosSB,mulCosBB;
    tmp = calloc(numOfData*numOfFeature,sizeof(double));
    for(i=0;i<numOfFeature/2;i++){
        for(j=0;j<numOfFeature/2;j++){
            tmp[i][j]=eigenvalue[i][j];
            theta=1/2*atan2(2*eigenvalue[i][j],eigenvalue[j][j]-eigenvalue[i][i]);
        }
        valueOfCos=cos(theta);
        valueOfSin=sin(theta);
        mulCosSS=valueOfCos*eigenvalue[small][small];
        mulCosSB=valueOfCos*eigenvalue[small][big];
        mulCosBB=valueOfCos*eigenvalue[big][big];
        mulSinSS=valueOfSin*eigenvalue[small][small];
        mulSinSB=valueOfSin*eigenvalue[small][big];
        mulSinBB=valueOfSin*eigenvalue[big][big];
        tmp[small][small]=valueOfCos*(mulCosSS-mulSinSB)-valueOfSin*(mulCosSB-mulSinBB);
        tmp[small][big]=valueOfSin*(mulCosSS-mulSinSB)+valueOfCos*(mulCosSB-mulSinBB);
        tmp[big][small]=tmp[small][big];
        tmp[big][big]=valueOfSin*(mulSinSS+mulCosSB)+valueOfCos*(mulSinSB+mulCosBB);
    }
    for(i=0;i<numOfFeature/2;i++){
        for(j=0;j<numOfFeature/2;j++){
            eigenvalue[i][j]=tmp[i][j];
            if(i==j){
            }else{
                
            }
        }
    }
    //計算終了：１　計算継続：０をreturn
}

//課題4 分母の固有値に+b（定数）を足す必要がある
