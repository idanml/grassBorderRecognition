#include "clustsigma.h"

clustSigma::clustSigma()
{
    this->sigma.resize(3,3);
    this->m=3;
    this->beta=0.5;
}


void clustSigma::sigmaCalc(matrix<double> var) {

    // Sigma Creationg
    double temp;
    temp=m*tgamma(m/(2*beta))/(pow(2,1/beta)*tgamma((m+2)/(2*beta)));
    this->sigma=temp*var;
    // End of Sigma Creation
}

matrix<double> clustSigma::getSigma(){
    return this->sigma;
}

