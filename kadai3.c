//4J02 s15015 �r�����i
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define numOfFeature 196        //������
#define numOfData 180           //�����f�[�^��
#define zero 0.01               //���R�r�@[���x]
#define maxNumberOfCalc 100   //���R�r�@[�ő�v�Z��]
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
    
    //�ۑ�3�ł́Csigmaxx.txt����f�[�^��ǂݎ�藘�p���邱�ƂƂ���
    i=0;
//    for(i=0;i<46;i++){
        sprintf(fnRead,"./sigmaData/%s%02d.txt",fnReadFormat,i+1);
        readData(eigenvalue,fnRead,numOfFeature,numOfFeature);    //eigenvalue ��sigma�f�[�^������
        sprintf(fnWrite,"./vectorData/%s%02d.txt",fnWriteFormat,i+1);
        calcEigenvalue(eigenvalue,eigenvector);
        writeDataTwoDim(eigenvalue,fnWrite);
//    }
}

//�t�@�C������f�[�^��ǂݎ��(�s���C�񐔎w��\)
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

//�s���\������(�z��́C196*196)
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

//�t�@�C���ɔz����������ށi�z��́C180�j
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

//�t�@�C���ɂQ�����z�����������(�z��́C196*196)
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
double oneElement_calcFeature(double data[][numOfFeature],int target)
{
    int i;
    double result=0;
    for(i=0;i<numOfData;i++){
        result+=data[i][target];
    }
    return result/numOfData;
}

//�w�萔�̕������̕��ϓ����ʂ��v�Z(���)
void average_calcFeature(double data[][numOfFeature],double average[])
{
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

//�ŗL�l�E�ŗL�x�N�g���̌v�Z
//(196����+1)�̓_�����ʂ��P���܂�
//�������A�����1?180�̃f�[�^�̂ݗ��p����
//�{���ł���΁A197�Ȃ��ƌŗL�l�E�ŗL�x�N�g���͋��܂�Ȃ����A�v�Z�덷�̂��߁A���܂�
//�덷������l�Ȃ̂ŁA���K�����n��e1~e180�܂ł𗘗p���邱�Ƃɂ���i��{�I�ɂ́j
//�����A���܂�Ȃ��ꍇ�����邽�߁A���܂�Ȃ��Ȃ������_�Ōv�Z��ł��؂�K�v������

//�ŗL�l�̌v�Z�����ɂ��āA�v�Z��ł��؂�Ȃǂ̔�����s��
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
        matrixP[i][i]=1;    //�Ίp�s��̍쐬
        eigenvector[i][i]=1;
    }
    start = clock();
    while(calcEigenvalueExe(eigenvalue,eigenvector,matrixP,matrixTmp)==1){
        end = clock();
        printf("  ���s���ԁF%.2f�b\n",(double)(end-start)/CLOCKS_PER_SEC);
    }
    //dispMatrix(eigenvalue);
    free(matrixP);
    free(matrixTmp);
}

//�ŗL�l�̌v�Z������1��s��
int calcEigenvalueExe(double eigenvalue[][numOfFeature],double eigenvector[][numOfFeature],double matrixP[][numOfFeature],double matrixTmp[][numOfFeature])
{
    int i,j;
    int small,big;  //(i,j)��(small,big)�ŕ\��
    /*�s��P
     *      ____small____ big ____
     *      ______________________
     * small____ cos��____ sin��____
     *      ______________________
     *  big ____-sin��____ cos��____
     *      ______________________
     * small<big
    */
    static int counter=0;
    position_k position;
    counter++;
    
    //�s��A(�ŏI�I�ɂ͌ŗL�l�ƂȂ�j�Feigenvalue�̔�Ίp�����̂�����ԑ傫�Ȓl��T��
    position = serchBiggest(eigenvalue);
    printf("small=%d big=%d[%.4lf]",position.row,position.column,eigenvalue[position.row][position.column]);

    //�֐�serchBiggest�ł́C�E��̕�����T���̂ŁCsmall,big�͂��̂悤�ɂȂ�
    //small => position.row  big => position.column
    createMatrix(matrixP,eigenvalue,position.row,position.column);

    //�ŗL�x�N�g�����X�V���� �ŗL�x�N�g���́C��ԑ傫���ŗL�l�̂��̂ɂȂ�
    //��P_1 * P_2 * P_3 .... P_n-1 * P_n
    calcProduct(eigenvector,matrixP,matrixTmp);
    copyMatrix(matrixTmp,eigenvector);
    
    //�ŗL�l���X�V����
    //��P^-1 * A * P = �����v�Z���� �ŗL�l�F��
    //A*P
    calcProduct(eigenvalue,matrixP,matrixTmp);
    copyMatrix(matrixTmp,eigenvalue);
    //P^-1 * A    
    turnMatrix(matrixP,position.row,position.column);
    calcProduct(matrixP,eigenvalue,matrixTmp);
    copyMatrix(matrixTmp,eigenvalue);

    //�v�Z�I���F-1�@�v�Z�p���F1��return
    //****************************************
    //���̉���if�������ɁC��Ίp������0�ɂȂ����Ƃ��̐���I�������������
    //****************************************
    if(counter>=maxNumberOfCalc){    //�v�Z�����K��l�𒴂����ꍇ�C�v�Z�I��
        printf("�v�Z���~\n");
        return -1;
    }else{
        printf("�v�Z��:%d",counter);
        return 1;
    }
}

//�s�񓯎m�̊|���Z���s��
//�s��left * �s��right 
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

//��Ίp�����̂����C�ő�l�̍s�ԍ��C��ԍ���Ԃ��֐�
//�Ώۍs��̂��߁C�s��̉E��̂ݒT�����邱�ƂƂ���
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
            if(max<fabs(target[i][j])){ //��Βl�Ŕ�r����
                max=target[i][j];
                position.row=i;
                position.column=j;
            }
        }
    }
    return position;    //��Ίp�����̂����C�ő�l�̍s�ԍ��C��ԍ���Ԃ�
}

//�s��P���쐬����
//�s��P�́C�ŏ��ɕʊ֐����i�֐�oneWord_calcEigenvalue()�j�ɂĐ錾����matrixP�𗘗p����
//matrixP�́C�P�ʍs��ł���C������X�V���Ďg���܂킷���ƂƂ���
//���̂��߁C�O��̏ꏊ��static�ɕێ����Ēu���C�V�����l�ɍX�V����ۂɁC������x�P�ʍs��ɖ߂��������s��
void createMatrix(double matrixP[][numOfFeature],double eigenvalue[][numOfFeature],int small,int big)
{
    static int pastSmall=0,pastBig=1;    //�O��̈ʒu��ۑ��B����(0,1)�Ȃ̂́C(0,0)�ɂ����1�s1��̒l���������邽��
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

//�s��P��]�u���C�s��P^-1�ɂ��鏈��
//�s��P�͑Ώۍs��ł���C�ω�������̂́Csin��-sin���������邾���ł悢
void turnMatrix(double target[][numOfFeature],int small,int big)
{ 
    target[small][big]=target[small][big]*(-1.0);
    target[big][small]=target[big][small]*(-1.0);
}

//�z��̃R�s�[
void copyMatrix(double origin[][numOfFeature],double target[][numOfFeature])
{
    int i,j;
    for(i=0;i<numOfFeature;i++){
        for(j=0;j<numOfFeature;j++){
            target[i][j]=origin[i][j];
        }
    }
}

//�Ίp�s�񂩂ǂ����m�F����
//�s��target���Ώۍs��ł���O��Ȃ̂ŁC�E��̕��������ׂă[���ł��邩�m�F���s����
//���[������́C���炩���ߒ�`����zero�����ł��邩�Ŋm�F
//�Ԃ�l�@�Ίp�s��F1�@��Ίp�s��F-1
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