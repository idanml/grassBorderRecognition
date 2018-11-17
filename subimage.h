#ifndef SUBIMAGE_H
#define SUBIMAGE_H
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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost::numeric::ublas;

class subImage {
    QImage sImg;
    matrix<double> red;
    matrix<double> green;
    matrix<double> blue;
    matrix<double> waveRed;
    matrix<double> waveGreen;
    matrix<double> waveBlue;
    matrix<double> normRed;
    matrix<double> normGreen;
    matrix<double> normBlue;
    matrix<double> sigma;
    double m, beta;

  public:
    subImage ();
    setData(int h, int w, QImage img);
    QImage getImg();
    void setRed(double num ,int i,int j);
    void setGreen(double num,int i,int j);
    void setBlue(double num ,int i,int j);
    void matrixToArray (double* redArr , double* greenArr , double* blueArr );
    void matrixToNormArray (double* redArr , double* greenArr , double* blueArr );
    void matrixToWaveArray (double* redArr , double* greenArr , double* blueArr );
    void arrayToWaveMatrix (double redArr[] , double greenArr[] , double blueArr[]);
    void arrayToNormMatrix (double redArr[] , double greenArr[] , double blueArr[]);
    double getM();
    double getBeta();
    matrix<double> getSigma();
    void setSigma(matrix<double> sigma);
    double getWaveRed(int i, int j);
    double getWaveGreen(int i, int j);
    double getWaveBlue(int i, int j);
    double getRed(int i, int j);
    double getGreen(int i, int j);
    double getBlue(int i, int j);
    double getNormRed(int i, int j);
    double getNormGreen(int i, int j);
    double getNormBlue(int i, int j);
    void testImg (QPoint index, QColor color);
};
#endif // SUBIMAGE_H
