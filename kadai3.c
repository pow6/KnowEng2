//4J02 s15015 池口恭司
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define numOfFeature 196        //特徴量
#define numOfData 180           //文字データ数
#define zero 0.01               //ヤコビ法[精度]
#define maxNumberOfCalc 100   //ヤコビ法[最大計算量]
#define fileRead(fileName,fileStream) fileStream=fopen(fileName,"r");if(fileStream==NULL){printf("cannat read file[%s]",fileName);exit(1);}
#define fileWrite(fileName,fileStream) fileStream=fopen(fileName,"w");if(fileStream==NULL){printf("cannat write file[%s]",fileName);exit(1);}

typedef struct {
    int row;
    int column;
} position_k;

void readData(double data[][numOfFeature],char fileName[],int row,int column);
void dispMatrix(double data[][numOfFeature]);
void writeData(double average[],char fileName[]);
void writeDataTwoDim(double coriance[][numOfFeature],char fileName[]);
double oneElement_calcFeature(double data[][numOfFeature],int target);
void average_calcFeature(double data[][numOfFeature],double average[]);
void calcCovariance(double data[][numOfFeature],double average[],double covariance[][numOfFeature]);
void calcEigenvalue(double eigenvalue[][numOfFeature],double eigenvector[][numOfFeature]);
int calcEigenvalueExe(double eigenvalue[][numOfFeature],double eigenvector[][numOfFeature],double matrixP[][numOfFeature],double matrixTmp[][numOfFeature]);
void calcProduct(double left[][numOfFeature],double right[][numOfFeature],double result[][numOfFeature]);
position_k serchBiggest(double target[][numOfFeature]);
void createMatrix(double matrixP[][numOfFeature],double eigenvalue[][numOfFeature],int small,int big);
void turnMatrix(double target[][numOfFeature],int small,int big);
void copyMatrix(double origin[][numOfFeature],double target[][numOfFeature]);
int judgeDiagonal(double target[][numOfFeature]);

void main()
{
    char fnReadFormat[12]="sigma",fnWriteFormat[12]="vector";
    char fnRead[24],fnWrite[24];
    int i,j;
    double (*eigenvalue)[numOfFeature],(*eigenvector)[numOfFeature];
    eigenvalue = malloc(sizeof(double)*numOfFeature*numOfFeature);
    eigenvector = malloc(sizeof(double)*numOfFeature*numOfFeature);
    
    //課題3では，sigmaxx.txtからデータを読み取り利用することとする
    i=0;
//    for(i=0;i<46;i++){
        sprintf(fnRead,"./sigmaData/%s%02d.txt",fnReadFormat,i+1);
        readData(eigenvalue,fnRead,numOfFeature,numOfFeature);    //eigenvalue にsigmaデータを入れる
        sprintf(fnWrite,"./vectorData/%s%02d.txt",fnWriteFormat,i+1);
        calcEigenvalue(eigenvalue,eigenvector);
        writeDataTwoDim(eigenvalue,fnWrite);
//    }
}

//ファイルからデータを読み取る(行数，列数指定可能)
void readData(double data[][numOfFeature],char fileName[],int row,int column)
{
    FILE *read;
    int i,j;
    fileRead(fileName,read);
    printf("Read[%s]\n",fileName);
    for(j=0;j<row;j++){
        for(i=0;i<column;i++){
            fscanf(read,"%lf",&data[j][i]);
        }
    }
    fclose(read);
}

//行列を表示する(配列は，196*196)
void dispMatrix(double data[][numOfFeature])
{
    int i,j;
    for(i=0;i<numOfFeature;i++){
        for(j=0;j<numOfFeature;j++){
            printf("%.8lf ,",data[i][j]);
        }
        printf("\n");
    }
}

//ファイルに配列を書き込む（配列は，180）
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

//ファイルに２次元配列を書き込む(配列は，196*196)
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

//固有値・固有ベクトルの計算
//(196次元+1)個の点→平面が１つ決まる
//しかし、今回は1?180個のデータのみ利用する
//本来であれば、197個ないと固有値・固有ベクトルは求まらないが、計算誤差のため、求まる
//誤差がある値なので、正規直交系はe1~e180までを利用することにする（基本的には）
//ただ、求まらない場合もあるため、求まらなくなった時点で計算を打ち切る必要がある

//固有値の計算処理について、計算を打ち切るなどの判定を行う
void calcEigenvalue(double eigenvalue[][numOfFeature],double eigenvector[][numOfFeature])
{
    clock_t start,end;
    int i,j;
    double (*matrixP)[numOfFeature],(*matrixTmp)[numOfFeature];
    matrixP = calloc(numOfFeature*numOfFeature,sizeof(double));
    matrixTmp = malloc(sizeof(double)*numOfFeature*numOfFeature);
    for(i=0;i<numOfFeature;i++){
        for(j=0;j<numOfFeature;j++){
            eigenvector[i][j]=0;
        }
        matrixP[i][i]=1;    //対角行列の作成
        eigenvector[i][i]=1;
    }
    start = clock();
    while(calcEigenvalueExe(eigenvalue,eigenvector,matrixP,matrixTmp)==1){
        end = clock();
        printf("  実行時間：%.2f秒\n",(double)(end-start)/CLOCKS_PER_SEC);
    }
    //dispMatrix(eigenvalue);
    free(matrixP);
    free(matrixTmp);
}

//固有値の計算処理を1回行う
int calcEigenvalueExe(double eigenvalue[][numOfFeature],double eigenvector[][numOfFeature],double matrixP[][numOfFeature],double matrixTmp[][numOfFeature])
{
    int i,j;
    int small,big;  //(i,j)を(small,big)で表す
    /*行列P
     *      ____small____ big ____
     *      ______________________
     * small____ cosθ____ sinΘ____
     *      ______________________
     *  big ____-sinθ____ cosΘ____
     *      ______________________
     * small<big
    */
    static int counter=0;
    position_k position;
    counter++;
    
    //行列A(最終的には固有値となる）：eigenvalueの非対角成分のうち一番大きな値を探す
    position = serchBiggest(eigenvalue);
    printf("small=%d big=%d[%.4lf]",position.row,position.column,eigenvalue[position.row][position.column]);

    //関数serchBiggestでは，右上の部分を探すので，small,bigはこのようになる
    //small => position.row  big => position.column
    createMatrix(matrixP,eigenvalue,position.row,position.column);

    //固有ベクトルを更新する 固有ベクトルは，一番大きい固有値のものになる
    //⇒P_1 * P_2 * P_3 .... P_n-1 * P_n
    calcProduct(eigenvector,matrixP,matrixTmp);
    copyMatrix(matrixTmp,eigenvector);
    
    //固有値を更新する
    //⇒P^-1 * A * P = Λを計算する 固有値：Λ
    //A*P
    calcProduct(eigenvalue,matrixP,matrixTmp);
    copyMatrix(matrixTmp,eigenvalue);
    //P^-1 * A    
    turnMatrix(matrixP,position.row,position.column);
    calcProduct(matrixP,eigenvalue,matrixTmp);
    copyMatrix(matrixTmp,eigenvalue);

    //計算終了：-1　計算継続：1をreturn
    //****************************************
    //この下のif文処理に，非対角成分が0になったときの正常終了処理も入れる
    //****************************************
    if(counter>=maxNumberOfCalc){    //計算数が規定値を超えた場合，計算終了
        printf("計算中止\n");
        return -1;
    }else{
        printf("計算数:%d",counter);
        return 1;
    }
}

//行列同士の掛け算を行う
//行列left * 行列right 
void calcProduct(double left[][numOfFeature],double right[][numOfFeature],double result[][numOfFeature])   
{
    int i,j,z;
    for(i=0;i<numOfFeature;i++){
        for(j=0;j<numOfFeature;j++){
            result[i][j]=0;
            for(z=0;z<numOfFeature;z++){
                result[i][j]+=left[i][z]*right[z][j];
            }
        }
    }
}

//非対角成分のうち，最大値の行番号，列番号を返す関数
//対象行列のため，行列の右上のみ探索することとする
position_k serchBiggest(double target[][numOfFeature])
{
    position_k position;
    int i,j;
    double max;
    max=target[0][1];
    position.row=0;
    position.column=1;
    for(i=0;i<numOfFeature-1;i++){
        for(j=i+1;j<numOfFeature;j++){
            if(max<fabs(target[i][j])){ //絶対値で比較する
                max=target[i][j];
                position.row=i;
                position.column=j;
            }
        }
    }
    return position;    //非対角成分のうち，最大値の行番号，列番号を返す
}

//行列Pを作成する
//行列Pは，最初に別関数内（関数oneWord_calcEigenvalue()）にて宣言したmatrixPを利用する
//matrixPは，単位行列であり，それを更新して使いまわすこととする
//そのため，前回の場所をstaticに保持して置き，新しい値に更新する際に，もう一度単位行列に戻す処理を行う
void createMatrix(double matrixP[][numOfFeature],double eigenvalue[][numOfFeature],int small,int big)
{
    static int pastSmall=0,pastBig=1;    //前回の位置を保存。初回(0,1)なのは，(0,0)にすると1行1列の値が消失するため
    double theta;
    matrixP[pastSmall][pastSmall]=1;
    matrixP[pastBig][pastBig]=1;
    matrixP[pastSmall][pastBig]=0;
    matrixP[pastBig][pastSmall]=0;
    theta=0.5*atan(2*eigenvalue[small][big]/(eigenvalue[big][big]-eigenvalue[small][small]));
    printf(" theta=%.4lf ass=%.4lf abb=%.4lf abs=%.4lf asb=%.4lf ",theta,eigenvalue[small][small],eigenvalue[big][big],eigenvalue[big][small],eigenvalue[small][big]);
    matrixP[small][small]=cos(theta);
    matrixP[big][big]=matrixP[small][small];
    matrixP[small][big]=sin(theta);
    matrixP[big][small]=matrixP[small][big]*(-1);
    pastSmall = small;
    pastBig = big;
}

//行列Pを転置し，行列P^-1にする処理
//行列Pは対象行列であり，変化させるのは，sinと-sinを交換するだけでよい
void turnMatrix(double target[][numOfFeature],int small,int big)
{ 
    target[small][big]=target[small][big]*(-1.0);
    target[big][small]=target[big][small]*(-1.0);
}

//配列のコピー
void copyMatrix(double origin[][numOfFeature],double target[][numOfFeature])
{
    int i,j;
    for(i=0;i<numOfFeature;i++){
        for(j=0;j<numOfFeature;j++){
            target[i][j]=origin[i][j];
        }
    }
}

//対角行列かどうか確認する
//行列targetが対象行列である前提なので，右上の部分がすべてゼロであるか確認を行う※
//※ゼロ判定は，あらかじめ定義したzero未満であるかで確認
//返り値　対角行列：1　非対角行列：-1
int judgeDiagonal(double target[][numOfFeature])
{
    int i,j;
    for(i=0;i<numOfFeature;i++){
        for(j=i+1;j<numOfFeature;j++){
            if(target[i][j]>zero){
                return -1;
            }
        }
    }
    return 1;
}