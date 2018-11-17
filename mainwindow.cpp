#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QCoreApplication>
#include <QRect>
#include <QtGui/QImage>
#include <QApplication>
#include <QtWidgets>
#include <iostream>
#include <sstream>
#include <QVector>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "haar.hpp"
#include "subimage.h"
#include "mainwindow.h"
#include "clustsigma.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <chrono>
#include <ctime>

extern "C" {
#include <src/cluster.h>
}

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
using namespace std;
using namespace boost::numeric::ublas;
using Eigen::MatrixXd;
using Eigen::VectorXcd;

//History Node
struct node {
    QString time;
    QString name;
    QString persent;
    int feedback;
    node* next;
};

void initNode(node *head){
    head->next = NULL;
}

void insertFront(node **head, QString time, QString fileName, QString persent, int feedback) {
    node *newNode = new node;
    newNode->time = time;
    newNode->name = fileName;
    newNode->persent = persent;
    newNode->next = *head;
    *head = newNode;
}


int calcSucc (node *head) {
    // Calaculate the success rate of the program
    node *cur = head;
    double count = 0, succ = 0;

    while(cur->next) {
        count ++;
        if (cur->feedback == 1)
            succ++;

        cur = cur->next;
    }

    double res = succ/count;
    return (((res)*100));

}

void normalize(double waveRed[],double waveGreen[],double waveBlue[],double* normRed, double* normGreen,double* normBlue,int arraySize, double level)
{
    double dev=0, avg=0;
    //Deviation calculation
             // Calc Avg
             for (int p=0 ; p<arraySize ; p++ ) {
                avg+=(waveRed[p]+waveGreen[p]+waveBlue[p]);
             }
             avg=avg/(arraySize*3);

             // End of AVG calc


             for (int p=0 ; p<arraySize ; p++ ) {
                 dev+=(pow((waveRed[p]-avg),2)+pow((waveGreen[p]-avg),2)+pow((waveBlue[p]-avg),2));
             }
             dev=sqrt(dev/(arraySize*3));
    // End of Deviation calculation

}

void createMatrix (matrix<double> &var, matrix<subImage> subImg, int i, int j, int arraySize) {

    double rr=0, rg=0, rb=0, gg=0, bb=0, gb=0;

    double normRed[arraySize];
    double normGreen[arraySize];
    double normBlue[arraySize];

    subImg(i,j).matrixToWaveArray(normRed, normGreen, normBlue); //  // Change 2שים לב שאתה עובד פה עם המערכים שיצרת לפני שניה זה הכל זבל ולא ערכים שקיבלנו

    // Create Var Matrix
   for (int p=0 ; p<arraySize ; p++ ) {
       rr+=normRed[p]*normRed[p];
       rg+=normRed[p]*normGreen[p];
       rb+=normRed[p]*normBlue[p];
       gg+=normGreen[p]*normGreen[p];
       bb+=normBlue[p]*normBlue[p];
       gb+=normGreen[p]*normBlue[p];
   }

   var(0,0)=rr/arraySize;
   var(1,1)=gg/arraySize;
   var(2,2)=bb/arraySize;
   var(0,1)=var(1,0)=rg/arraySize;
   var(0,2)=var(2,0)=rb/arraySize;
   var(2,1)=var(1,2)=gb/arraySize;
    // End of Create Var Matrix
}

matrix<double> createSigma (matrix<subImage> subImg, int i, int j, matrix<double> var) {

    double temp, m, beta;
    matrix<double> sigma (3, 3);

    m=subImg(i,j).getM();
    beta=subImg(i,j).getBeta();

    // Sigma Creationg

    temp=m*tgamma(m/(2*beta))/(pow(2,1/beta)*tgamma((m+2)/(2*beta)));
    sigma=temp*var;
    return (sigma);
    // End of Sigma Creation
}

void distCalc(subImage subImgagesArr[],  matrix<double> &dist, double** distMat, int subPics)
{
    matrix<double> sigma1(3,3), sigma2(3,3), invSigma2(3,3), result(3,3);
    matrix<complex<double>> eigenVal(1,3), s1(1,1), s2(3,3);
    complex<double> sum1,sum2;
    double bh=4/(4*5);           //m=3

    for(int i=0; i<subPics; i++) {
        for(int j=i; j<subPics; j++)
        {
            if(i==j){
                dist(i,j)=0;
                distMat[i][j]=0;
            }
            else{
         sigma1=subImgagesArr[i].getSigma();
         sigma2=subImgagesArr[j].getSigma();
         double determinant =    sigma2(0,0)*(sigma2(1,1)*sigma2(2,2)-sigma2(2,1)*sigma2(1,2))
                                 -sigma2(0,1)*(sigma2(1,0)*sigma2(2,2)-sigma2(1,2)*sigma2(2,0))
                                 +sigma2(0,2)*(sigma2(1,0)*sigma2(2,1)-sigma2(1,1)*sigma2(2,0));
         double invdet = 1/determinant;
         invSigma2(0,0) =  (sigma2(1,1)*sigma2(2,2)-sigma2(2,1)*sigma2(1,2))*invdet;
         invSigma2(1,0) = -(sigma2(0,1)*sigma2(2,2)-sigma2(0,2)*sigma2(2,1))*invdet;
         invSigma2(2,0) =  (sigma2(0,1)*sigma2(1,2)-sigma2(0,2)*sigma2(1,1))*invdet;
         invSigma2(0,1) = -(sigma2(1,0)*sigma2(2,2)-sigma2(1,2)*sigma2(2,0))*invdet;
         invSigma2(1,1) =  (sigma2(0,0)*sigma2(2,2)-sigma2(0,2)*sigma2(2,0))*invdet;
         invSigma2(2,1) = -(sigma2(0,0)*sigma2(1,2)-sigma2(1,0)*sigma2(0,2))*invdet;
         invSigma2(0,2) =  (sigma2(1,0)*sigma2(2,1)-sigma2(2,0)*sigma2(1,1))*invdet;
         invSigma2(1,2) = -(sigma2(0,0)*sigma2(2,1)-sigma2(2,0)*sigma2(0,1))*invdet;
         invSigma2(2,2) =  (sigma2(0,0)*sigma2(1,1)-sigma2(1,0)*sigma2(0,1))*invdet;
         result=prod(sigma1,invSigma2);

         MatrixXd resulTemp = MatrixXd::Ones(3,3);

         for (int p = 0 ; p < 3 ; p++)
             for (int d = 0 ; d < 3 ; d++)
                 resulTemp(p,d)=result(p,d);
         VectorXcd eivals = resulTemp.eigenvalues();
         for (int t = 0 ; t < 3 ; t++)
             eigenVal(0,t)=log(eivals[t]);
         s1=prod(eigenVal,trans(eigenVal));
         sum1=s1(0,0);
         s2=prod(trans(eigenVal),eigenVal);
         sum2=(s2(0,0)+s2(0,1)+s2(0,2)+s2(1,0)+s2(1,1)+s2(1,2)+s2(2,0)+s2(2,1)+s2(2,2)-sum1)/2;
        dist(i,j)=dist(j,i)=distMat[i][j]=distMat[j][i]=pow(abs(((3*bh - 0.25) * sum1 + 2*(bh - 0.25) * sum2)),0.5);

            }
        }
    }
}

void graphCreation(Graph &g,int s,  Graph::vertex_descriptor v[], int clusterid[], int flag)
{
    for(int i=0; i<400; i++)
      v[i]=add_vertex(g);

    for(int i=0; i<380; i++)
    {
        if(s==clusterid[i] && clusterid[i]==clusterid[i+20]) //check the subpic below
            add_edge(v[i], v[i+20],g);
        if(s==clusterid[i] && clusterid[i]==clusterid[i+1] && i%20!=19) //check the subpic on the right side
            add_edge(v[i], v[i+1],g);
        if(flag==1)
          {
            if(s==clusterid[i] && clusterid[i]==clusterid[i+21] && i%20!=19) //check the right diag
                add_edge(v[i], v[i+21],g);
            if(s==clusterid[i] && clusterid[i]==clusterid[i+19] && i%20!=0) // cheack the left diag
                add_edge(v[i], v[i+19],g);
         }
    }
    for(int i=380; i<399; i++)
        if(s==clusterid[i] && clusterid[i]==clusterid[i+1]) //check the subpic on the right side
            add_edge(v[i], v[i+1],g);
}

int findMax(int num,std::vector<int> &component, int size, int *numOfSubpic)
 {
     int max, maxc;
     int count[num]={0};
     for(int i=0; i<size; i++)
        count[component[i]]++;
     max=count[0];
     maxc=0;
     for(int i=1; i<num; i++)
         if(max<count[i])
         {
             max= count[i];
             maxc=i;
         }
     *numOfSubpic=max;
     return maxc;
 }

void clustDistCalc(clustSigma clustArr[],  matrix<double> &dist, double** distMat, int numOfClust)
 {
     matrix<double> sigma1(3,3), sigma2(3,3), invSigma2(3,3), result(3,3);
     matrix<complex<double>> eigenVal(1,3), s1(1,1), s2(3,3);
     complex<double> sum1,sum2;
     double bh=4/(4*5);           //m=3

     for(int i=0; i<numOfClust; i++) {
         for(int j=i; j<numOfClust; j++)
         {
             if(i==j){
                 dist(i,j)=0;
                 distMat[i][j]=0;
             }
             else{
          sigma1=clustArr[i].getSigma();
          sigma2=clustArr[j].getSigma();
          double determinant =    sigma2(0,0)*(sigma2(1,1)*sigma2(2,2)-sigma2(2,1)*sigma2(1,2))
                                  -sigma2(0,1)*(sigma2(1,0)*sigma2(2,2)-sigma2(1,2)*sigma2(2,0))
                                  +sigma2(0,2)*(sigma2(1,0)*sigma2(2,1)-sigma2(1,1)*sigma2(2,0));
          double invdet = 1/determinant;
          invSigma2(0,0) =  (sigma2(1,1)*sigma2(2,2)-sigma2(2,1)*sigma2(1,2))*invdet;
          invSigma2(1,0) = -(sigma2(0,1)*sigma2(2,2)-sigma2(0,2)*sigma2(2,1))*invdet;
          invSigma2(2,0) =  (sigma2(0,1)*sigma2(1,2)-sigma2(0,2)*sigma2(1,1))*invdet;
          invSigma2(0,1) = -(sigma2(1,0)*sigma2(2,2)-sigma2(1,2)*sigma2(2,0))*invdet;
          invSigma2(1,1) =  (sigma2(0,0)*sigma2(2,2)-sigma2(0,2)*sigma2(2,0))*invdet;
          invSigma2(2,1) = -(sigma2(0,0)*sigma2(1,2)-sigma2(1,0)*sigma2(0,2))*invdet;
          invSigma2(0,2) =  (sigma2(1,0)*sigma2(2,1)-sigma2(2,0)*sigma2(1,1))*invdet;
          invSigma2(1,2) = -(sigma2(0,0)*sigma2(2,1)-sigma2(2,0)*sigma2(0,1))*invdet;
          invSigma2(2,2) =  (sigma2(0,0)*sigma2(1,1)-sigma2(1,0)*sigma2(0,1))*invdet;
          result=prod(sigma1,invSigma2);

          MatrixXd resulTemp = MatrixXd::Ones(3,3);

          for (int p = 0 ; p < 3 ; p++)
              for (int d = 0 ; d < 3 ; d++)
                  resulTemp(p,d)=result(p,d);
          VectorXcd eivals = resulTemp.eigenvalues();
          for (int t = 0 ; t < 3 ; t++)
              eigenVal(0,t)=log(eivals[t]);
          s1=prod(eigenVal,trans(eigenVal));
          sum1=s1(0,0);
          s2=prod(trans(eigenVal),eigenVal);
          sum2=(s2(0,0)+s2(0,1)+s2(0,2)+s2(1,0)+s2(1,1)+s2(1,2)+s2(2,0)+s2(2,1)+s2(2,2)-sum1)/2;
         dist(i,j)=dist(j,i)=distMat[i][j]=distMat[j][i]=pow(abs(((3*bh - 0.25) * sum1 + 2*(bh - 0.25) * sum2)),0.5);

             }
         }
     }
 }


void fixCluster (int *clustersArray, int clusterID, int nSubPic) {

     matrix<int> clustersIds (nSubPic, nSubPic);
     int c=0, check; // Counter

     // Array to Matrix
     for (int i=0 ; i < nSubPic ; i++)
         for (int j=0 ;j < nSubPic ; j++) {
             clustersIds(i,j)= clustersArray[c];
             c++;
         }
     // End of Array to Matrix

     // Fixing
     for (int i=0 ; i < clustersIds.size1() ; i++)
         for (int j=0 ;j < clustersIds.size2() ; j++) {
            if (clustersIds(i,j)!=clusterID)
            {
                check = 0;

                if (j!=0)
                {
                    if (clustersIds(i,j-1)==clusterID )
                        check++;
                }


                if (j!=clustersIds.size2()-1)
                {
                    if (clustersIds(i,j+1)==clusterID )
                        check++;
                }

                if (i!=0)
                {
                    if (clustersIds(i-1,j)==clusterID )
                        check++;
                }

                if (i!=clustersIds.size1()-1)
                {
                    if (clustersIds(i+1,j)==clusterID )
                        check++;
                }

                if (check>=3)
                {
                    clustersIds(i,j)=clusterID;
                    continue;
                }

                //Checking for corners

                // Top Left
                if (i==0 && j==0) {
                    if (clustersIds(i+1,j+1)==clusterID && clustersIds(i,j+1)==clusterID && clustersIds(i+1,j)==clusterID)
                        clustersIds(i,j)=clusterID;
                }

                //Top Right
                if (i==0 && j==clustersIds.size2()-1) {
                    if (clustersIds(i+1,j-1)==clusterID && clustersIds(i,j-1)==clusterID && clustersIds(i+1,j)==clusterID)
                        clustersIds(i,j)=clusterID;
                }

                //Bottom Left
                if (i==clustersIds.size1()-1 && j==0) {
                    if (clustersIds(i-1,j)==clusterID && clustersIds(i-1,j+1)==clusterID && clustersIds(i,j+1)==clusterID)
                        clustersIds(i,j)=clusterID;
                }

                //Bottom Right
                if (i==clustersIds.size1()-1 && j==clustersIds.size2()-1) {
                    if (clustersIds(i-1,j)==clusterID && clustersIds(i-1,j-1)==clusterID && clustersIds(i,j-1)==clusterID)
                        clustersIds(i,j)=clusterID;
                }
            }
         } // End of Fixing

     // Matrix to Array
     c=0;
     for (int i=0 ; i < clustersIds.size1() ; i++)
         for (int j=0 ;j < clustersIds.size2() ; j++) {
              clustersArray[c] = clustersIds(i,j);
             c++;
         }
     // End of Matrix to Array
}


node *head;
int flagInHistory = 0;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    head = new node;
    initNode(head);

    ui->setupUi(this);
    QStringList title;

    ui->stackedWidget->setCurrentIndex(0);
    ui->tableWidget->setColumnCount(4);
    title << "_____________Date_____________" << "_____________________________________Image_____________________________________" << "% of Grass" << "__Feedback__";
    ui -> tableWidget -> setHorizontalHeaderLabels(title);
    ui->tableWidget->resizeColumnsToContents();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_history_clicked()
{
    node *cur = head;
    int count = 0;


    if (flagInHistory==0) {
        ui->tableWidget->clearContents();
        for (int i = 0 ; i < ui->tableWidget->rowCount() ; i ++)
            ui->tableWidget->removeRow(i);

        while(cur->next) {
            ui->tableWidget->insertRow(count);
                //int count = ui->tableWidget->rowCount()-1;
                ui->tableWidget->setItem(count,DATE,new QTableWidgetItem(cur->time));
                ui->tableWidget->setItem(count,IMAGE,new QTableWidgetItem(cur->name));
                ui->tableWidget->setItem(count,PERSENT,new QTableWidgetItem(cur->persent));
                if (cur->feedback==1)
                    ui->tableWidget->setItem(count,FEEDBACK,new QTableWidgetItem("Correct"));
                if (cur->feedback==0)
                    ui->tableWidget->setItem(count,FEEDBACK,new QTableWidgetItem("Incorrect"));
                count++;
            cur = cur->next;
        }
    }


    flagInHistory = 1;
    ui->stackedWidget->setCurrentIndex(1);
}

void MainWindow::on_pushButton_main_clicked()
{
    flagInHistory = 0;



    // Calculate success rate
    if (head->next!=NULL) {
        int sRate = calcSucc(head);
        QString start = "Success Rate: ", end = "%";
        QString updatedRate = QString::fromStdString(boost::lexical_cast<std::string>(sRate));

        start.append(updatedRate);
        start.append(end);

        ui->successRate_label->setText(start);
    }
    // End of Calculate success rate

    ui->stackedWidget->setCurrentIndex(0);

}

void MainWindow::on_analize_btn_clicked(){
    flagInHistory = 0;

    QDateTime now = QDateTime::currentDateTime();
    QString time = now.toString();

   ui->browse_btn->setText("Browse");

    QImage img(fileName);  // Load an image

    // Algorithem

    // Declerations
    int nSubPic = 20, subPics = nSubPic*nSubPic;
    int img_w = img.width(), w = img_w/nSubPic;  //Sizes of image and subimages
    int img_h = img.height(), h = img_h/nSubPic;

    matrix<subImage> subImgages (nSubPic, nSubPic);       // Matrix for sub Images
    subImage subImgagesArr [subPics];
    int count=0;
    QRect rect(0, 0, w, h);
    QImage cropped;
    double max, min, level=2, dispRed[w*h], dispGreen[w*h], dispBlue[w*h],
            tempRed[w*h], tempGreen[w*h], tempBlue[w*h], normRed[w*h], normGreen[w*h], normBlue[w*h];
    // End of Declerations


    // Subimages Loop
    for(unsigned i = 0, k = 0; i < subImgages.size1(); i++){
        for(unsigned j = 0; j < subImgages.size2(); j++){

            rect.moveTo(w*j, h*i);
            cropped = img.copy(rect);  //exstracting subimages


//            name.append(QString::number(k)); // name of subimage
//            name.append(".jpg");

//            cropped.save(name,"JPG"); //saving subimages
//            k++;

            subImgages(i,j).setData(h, w, cropped);

            // RGB Creation
            for (int z = 0 ; z < h; ++z){
                for (int f = 0; f < w; ++f) {
                     QRgb c = cropped.pixel(f,z);
                     subImgages(i,j).setRed(qRed(c),z,f);
                     subImgages(i,j).setGreen(qGreen(c),z,f);
                     subImgages(i,j).setBlue(qBlue(c),z,f);
                   }
            }
            // End of RGB Creation

            subImgages(i,j).matrixToArray(tempRed, tempGreen, tempBlue);

            haar_2d (w,h,tempRed);
            haar_2d (w,h,tempGreen);
            haar_2d (w,h,tempBlue);

            subImgages(i,j).arrayToWaveMatrix(tempRed, tempGreen, tempBlue);


            //##################################################################
            int arraySize=w*h;       // Size of the arrays
            //normalize(tempRed,tempGreen,tempBlue,normRed,normGreen,normBlue,arraySize,level);

            subImgages(i,j).arrayToNormMatrix(normRed, normGreen, normBlue);

            matrix<double> var (3, 3);
            createMatrix (var, subImgages , i, j, arraySize);

            subImgages(i,j).setSigma(createSigma(subImgages , i, j, var));
            // End of Gauss

            subImgagesArr[count]=subImgages(i,j);
            count++;
        }
    }

    //End Subimages Loop

    //************************************************
    cout <<subImgagesArr[2].getSigma() << endl;  //[3,3]((273.946,539.369,797.244),(539.369,1062.26,1570.31),(797.244,1570.31,2321.48))
    cout << subImgagesArr[150].getSigma() << endl; //[3,3]((3393.92,3413.71,3438.12),(3413.71,3438.16,3466.4),(3438.12,3466.4,3499.7))
    //************************************************

    matrix<double> distances (subPics, subPics);

    // distMat initial
    double **distMat;
    distMat = new double *[subPics];

    for(int i = 0; i <subPics; i++)
        distMat[i] = new double[subPics];
    // END of initial

    distCalc(subImgagesArr, distances, distMat, subPics);
    //******************************************************
    cout <<  distMat[2][10] << endl; // 4.62053
    cout <<  distMat[15][19] << endl; // 2.53985
    cout <<  distMat[130][130] << endl; // 0
    //*****************************************************

    int clusterid[subPics], clusteridCopy[subPics];

    double *error;
    int *ifound;
    QColor col;

    error = new double;
    ifound = new int;
    int numOfClust = 4, k=0, clusters[numOfClust] = {0};

    kmedoids (numOfClust, subPics, distMat, 10, clusterid, error, ifound); // Clustering
    memcpy (clusteridCopy, clusterid, subPics*sizeof(int));
    sort(clusteridCopy,clusteridCopy+subPics);

    clusters[0]=clusteridCopy[0];

    for (int i=1; i < subPics ; i++) {
        if (clusteridCopy[i]!=clusters[k]) {
            k++;
            clusters[k]=clusteridCopy[i];

        }
    }
    fixCluster(clusterid,clusters[3],nSubPic);

    //**************************************************
    cout << clusters[0] << "  " << clusters[1] << "  " <<clusters[2] << "  " <<clusters[3] <<endl; //41  76  126  209
    //***************************************************

    int flag, numOfComp, maxComp, numOfSubpic, b;
    flag=1;     // if flag=1 check the subpic on the diagonal
    std::vector<int> component(subPics);
    Graph g[numOfClust];
    Graph::vertex_descriptor v[subPics];
    matrix<double> tempMat (3, 1);
    matrix<double> variMat (3, 3);
    clustSigma clSigma[numOfClust];



    for(int j=0;j<numOfClust; j++)
      {
          b=0;
          graphCreation(g[j], clusters[j] ,v ,clusterid ,flag);
          numOfComp = connected_components(g[j], &component[0]);
          maxComp=findMax(numOfComp,component, subPics, &numOfSubpic);
          //******************************************************
          cout << j << "  " << numOfComp << "  " << maxComp << "  " << numOfSubpic << endl; // 0  334  0  51
                                                                                            // 1  312  157  43
                                                                                            // 2  337  120  64
                                                                                            // 3  232  169  155
          // *****************************************************
          tempMat.resize(3, w*h*numOfSubpic, true);
          for(int i=0; i<subPics; i++)
          {
              if(maxComp==component[i]){
                  for(int l=0; l<w*h; l++)
                  {
                     tempMat(0,b*w*h+l)=subImgagesArr[i].getWaveRed(i/20, i%20);
                     tempMat(1,b*w*h+l)=subImgagesArr[i].getWaveGreen(i/20, i%20);
                     tempMat(2,b*w*h+l)=subImgagesArr[i].getWaveBlue(i/20, i%20);
              }
              b++;
              }
          }
          variMat=prod(tempMat, trans(tempMat))/(w*h*numOfSubpic);
          clSigma[j].sigmaCalc(variMat);
          //***************************************************
          cout << j << "   " << clSigma[j].getSigma() << endl; // 0   [3,3]((2258.16,4627.75,7148.46),(4627.75,9488.71,14657.6),(7148.46,14657.6,22644.2))
                                                               // 1   [3,3]((24.6121,26.8156,17.8537),(26.8156,29.745,19.2088),(17.8537,19.2088,14.0285))
                                                               // 2   [3,3]((44.2404,33.2085,20.7327),(33.2085,27.1092,16.8184),(20.7327,16.8184,14.1915))
                                                               // 3   [3,3]((539.821,857.474,305.22),(857.474,1372.87,479.941),(305.22,479.941,181.106))
          //***************************************************
      }

    matrix<double> clustDist (numOfClust, numOfClust);
    double **distClust;
    distClust = new double *[numOfClust];
    for(int i = 0; i <numOfClust; i++)
        distClust[i] = new double[numOfClust];
    clustDistCalc(clSigma, clustDist, distClust, numOfClust);
    //***********************************************************
    cout << clustDist << endl;  //[4,4]((0,3.75969,2.07886,0.234064),(3.75969,0,1.68083,3.99376),(2.07886,1.68083,0,2.31292),(0.234064,3.99376,2.31292,0))
    //**************************************************

    double score[2]; //numOfClust-2
    int clusteri[16]; //numOfClust*numOfClust
    for(int i=0 ;i<2; i++)
    {
       kmedoids (i+2, numOfClust, distClust, 10, clusteri, error, ifound);
       score[i]=*error;
       cout << score[i] << endl;
    }



    for(unsigned i = 0; i < subImgages.size1(); i++) {
            for(unsigned j = 0; j < subImgages.size2(); j++)
            {
                int mark = clusterid[i*subImgages.size1()+j];

                for (int z = 0 ; z < cropped.height(); z++) {
                    for (int f = 0; f < cropped.width(); f++) {
                         QPoint p(f,z);


                         if (mark==clusters[0])
                             col.setRgb(0, subImgages(i,j).getGreen(z,f), 0);
                         if (mark==clusters[1])
                             col.setRgb(subImgages(i,j).getRed(z,f), 0,0);
                         if (mark==clusters[2])
                             col.setRgb(subImgages(i,j).getRed(z,f), 0, 0);
                         if (mark==clusters[3])
                             col.setRgb(subImgages(i,j).getRed(z,f), subImgages(i,j).getGreen(z,f), subImgages(i,j).getBlue(z,f));

                         subImgages(i,j).testImg(p, col);
                    }
              }
            mark=-1;
    }
}


    QLabel *i_label[subPics];   // Labels on the window (places for subimages)

    int stratPos = 675-img_w/2;

    for(int i = 0, k = 0; i < nSubPic; i++){
        for(int j = 0; j < nSubPic; j++){
            rect.moveTo(w*j, h*i);
            QImage cropped = img.copy(rect);  //exstracting subimages

            i_label[k] = new QLabel("Final Photo", ui->newImage);  // add label to main window
            i_label[k]->setGeometry(stratPos+w*i, h*j, w, h); // set label position
            i_label[k]->setPixmap(QPixmap::fromImage(subImgages(j,i).getImg())); // load subimage to label
        }
    }
    insertFront(&head,time,fileName,fileName,-1);
    ui->stackedWidget->setCurrentIndex(2);
}




void MainWindow::on_browse_btn_clicked()
{
    fileName = QFileDialog::getOpenFileName(
                this, tr("Open File"), "C://", "Image files (*.png, *.jpg) ;; All files (*.*)"
                );
    ui->browse_btn->setText(fileName);


}

void MainWindow::on_correct_btn_clicked()
{
    head->feedback=1;
    on_pushButton_history_clicked();

}

void MainWindow::on_incorrect_btn_clicked()
{
    head->feedback=0;
    on_pushButton_history_clicked();

}
