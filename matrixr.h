#ifndef MATRIXR_H_
#define MATRIXR_H_

#include "vecr.h"

#define EIGEN_MATRIXBASE_PLUGIN "/home/eshamay/md/src/EigenMatrixAddon.h"
#include <Eigen/Core>
USING_PART_OF_NAMESPACE_EIGEN

#include <vector>

typedef enum {xx=0, yx=1, zx=2, xy=3, yy=4, zy=5, xz=6, yz=7, zz=8} element;

typedef Matrix3d	MatR;

#endif
