#include "subimage.h"

using namespace std;
using namespace boost::numeric::ublas;

subImage::subImage()
{
    this->m=3;
    this->beta=0.5;
    this->sigma.resize(m,m);
}

subImage::setData(int h, int w, QImage img)
{

    this->sImg=img;
    this->red.resize(h,w,true);
    this->green.resize(h,w,true);
    this->blue.resize(h,w,true);
    this->waveRed.resize(h,w,true);
    this->waveGreen.resize(h,w,true);
    this->waveBlue.resize(h,w,true);
    this->normRed.resize(h,w,true);
    this->normGreen.resize(h,w,true);
    this->normBlue.resize(h,w,true);
}


QImage subImage::getImg() {
    return (this->sImg);
}

void subImage::setRed(double num, int i, int j) {
    this->red(i,j)=num;
}

void subImage::setGreen(double num, int i, int j) {
    this->green(i,j)=num;
}

void subImage::setBlue(double num, int i, int j) {
    this->blue(i,j)=num;
}

double subImage::getRed(int i, int j) {
    return(this->red(i,j));
}

double subImage::getGreen(int i, int j) {
    return(this->green(i,j));
}

double subImage::getBlue(int i, int j) {
    return(this->blue(i,j));
}

double subImage::getWaveRed(int i, int j) {
    return(this->waveRed(i,j));
}

double subImage::getWaveGreen(int i, int j) {
    return(this->waveGreen(i,j));
}

double subImage::getWaveBlue(int i, int j) {
    return(this->waveBlue(i,j));
}

double subImage::getNormRed(int i, int j) {
    return(this->normRed(i,j));
}

double subImage::getNormGreen(int i, int j) {
    return(this->normGreen(i,j));
}

double subImage::getNormBlue(int i, int j) {
    return(this->normBlue(i,j));
}


double subImage::getM(){
    return this->m;
}

double subImage::getBeta(){
    return this->beta;
}

matrix<double> subImage::getSigma(){
    return this->sigma;
}

void subImage::setSigma(matrix<double> sigma){
    this->sigma=sigma;
}

void subImage::matrixToArray (double* redArr , double* greenArr , double* blueArr ) {
    int c = 0;

    for (int i=0 ; i<this->red.size1() ; i++) {
       for (int j=0 ; j<this->red.size2() ; j++) {
            redArr[c] = this->red(i,j);
            greenArr[c] = this->green(i,j);
            blueArr[c] = this->blue(i,j);
            c++;
       }
    }
}

void subImage::matrixToNormArray (double* redArr , double* greenArr , double* blueArr ) {
    int c = 0;

    for (int i=0 ; i<this->normRed.size1() ; i++) {
       for (int j=0 ; j<this->normRed.size2() ; j++) {
            redArr[c] = this->normRed(i,j);
            greenArr[c] = this->normGreen(i,j);
            blueArr[c] = this->normBlue(i,j);
            c++;
       }
    }
}

void subImage::matrixToWaveArray (double* redArr , double* greenArr , double* blueArr ) {
    int c = 0;

    for (int i=0 ; i<this->waveRed.size1() ; i++) {
       for (int j=0 ; j<this->waveRed.size2() ; j++) {
            redArr[c] = this->waveRed(i,j);
            greenArr[c] = this->waveGreen(i,j);
            blueArr[c] = this->waveBlue(i,j);
            c++;
       }
    }
}

void subImage::arrayToWaveMatrix (double redArr[] , double greenArr[] , double blueArr[] ) {
    int c = 0;

    for (int i=0 ; i<this->red.size1() ; i++) {
       for (int j=0 ; j<this->red.size2() ; j++) {
            this->waveRed(i,j)=redArr[c];
            this->waveGreen(i,j)=greenArr[c];
            this->waveBlue(i,j)=blueArr[c];
            c++;
       }
    }
}

void subImage::arrayToNormMatrix (double redArr[] , double greenArr[] , double blueArr[] ) {
    int c = 0;

    for (int i=0 ; i<this->normRed.size1() ; i++) {
       for (int j=0 ; j<this->normRed.size2() ; j++) {
            this->normRed(i,j)=redArr[c];
            this->normGreen(i,j)=greenArr[c];
            this->normBlue(i,j)=blueArr[c];
            c++;
       }
    }
}

void subImage::testImg (QPoint index, QColor color) {
    this->sImg.setPixelColor(index,color);
}
