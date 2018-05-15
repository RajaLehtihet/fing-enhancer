# Fingerprints image enhancer
Copyright (c)  2013-2015 Raja Lehtihet, All rights reserved.
## Introduction
This is a hybrid computational geometry - gray scale algorithm that estimates the ridge frequency of fingerprints image for enhancement. It generates a Delaunay Triangulation (DT) using extracted points of interest. This triangulation along with the local orientations give an accurate distance and orientation-based ridge frequency. A tuned anisotropic filter is locally applied and the enhanced output fingerprint image is obtained.

## Configuration and compilation
Using (qmake) the command line as follow: 
```
qmake
make -j4
```

## Usage
```
fing-enhance -f 101_2.tiff
```

## Algorithm
The proposed algorithm is composed of several stages:
- Block selection: Select the blocks that contain recoverable image parts. 
- Local minima extraction: Extract relevant points along the ridges. 
- Minima Delaunay Triangulation: Construct the DT using the extracted minima points. 
- Orientation estimation: Compute the dominant orientation of each block.
- Frequency estimation: Using the orientation and the DT, estimate the ridge frequency of each block. 
- Gabor filtering: Using the orientation, the ridge frequency of each block and its neighbouring blocks, we convolute the block with nine Gabor filters (9 being the eight-connected neighbouring blocks, plus the block itself). 
- Block equalization: Equalize the histogram of each final block.

## Prerequisites
* Qt5
* Libtiff
* Qhull

## Examples

![Original](https://github.com/RajaLehtihet/fing-enhancer/raw/master/images/orig.png)

The returned image presents salient ridges:

![Enhanced](https://github.com/RajaLehtihet/fing-enhancer/raw/master/images/enhanced.png)

## License
The present fingerprints image enhancer source code is licensed under Affero GPL version 3.0 or later. For any question  feel free to contact me at [Email](https://github.com/RajaLehtihet/fing-enhancer/raw/master/images/raja.png)

## Acknowledgements
Libtiff Tools Copyright (c) 2006, Rene Rebe, Copyright (c) 2004 for Andrey Kiselev
