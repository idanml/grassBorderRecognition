#ifndef CLUSTSIGMA_H
#define CLUSTSIGMA_H
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

class clustSigma
{
    matrix<double> sigma;
    double m, beta;
public:
    clustSigma();
    void sigmaCalc(matrix<double> var);
    matrix<double> getSigma();
};

#endif // CLUSTSIGMA_H
