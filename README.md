# grassBorderRecognition
Main program- mainwindow.cpp

Image processing program build for automatic lawn-mower, the main objective is to recognize grass in image by his texture.

The implemented algorithm:
Divide the image into smaller sub-images.

Wavelet coefficients is calculated for each sub-image.

Geodesic distance values between them evaluated.

Cluster the data and decide which sub-image contain grass.

Graph created and looking for connected components of sub-images (by there cluster).
