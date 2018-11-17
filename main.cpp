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
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>

extern "C" {
#include <src/cluster.h>
}



void rgbCreate (matrix<subImage> subImg, int i, int j) {
    // RGB Creation

    QImage cropped = subImg(i,j).getImg();

    for (int z = 0 ; z < cropped.height(); ++z) {
         int count=0;
           for (int f = 0; f < cropped.width(); ++f) {
             QRgb c = cropped.pixel(f,z);
             subImg(i,j).setRed(qRed(c),f,z);
             subImg(i,j).setGreen(qGreen(c),f,z);
             subImg(i,j).setBlue(qBlue(c),f,z);
             count++;
           }
     }
    // End of RGB Creation
}


int main(int argc, char *argv[]) {

//    QApplication app(argc, argv);  // Initializate an aplication

//    QWidget Main_Window;
//    Main_Window.resize(1024, 768);  // Create a main window

    QApplication app(argc, argv);
    MainWindow window;
    window.show();

// ###################################################################### Test our matrix ######################################################################
/*
    QLabel *i_label[subPics];   // Labels on the window (places for subimages)

    for(int i = 0, k = 0; i < nSubPic; i++){
        for(int j = 0; j < nSubPic; j++){
            rect.moveTo(w*j, h*i);
            QImage cropped = img.copy(rect);  //exstracting subimages

            i_label[k] = new QLabel("", &Main_Window);  // add label to main window
            i_label[k]->setGeometry(w*i, h*j, w, h); // set label position
            i_label[k]->setPixmap(QPixmap::fromImage(subImgages(j,i).getImg())); // load subimage to label
        }
    }
//###################################################################### Test our matrix ######################################################################
*/
//  Main_Window.show(); // Sow main windows
    return app.exec(); // return from application
}
