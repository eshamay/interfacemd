#pragma once
#ifndef VECR_H_
#define VECR_H_

#define EIGEN_MATRIXBASE_PLUGIN "/home/eshamay/md/src/EigenMatrixAddon.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
USING_PART_OF_NAMESPACE_EIGEN

#include <list>
 
typedef enum {x=0, y=1, z=2} coord;

typedef Eigen::Vector3d VecR;

typedef std::list< VecR, Eigen::aligned_allocator<VecR> > VecR_vec;
typedef VecR_vec::const_iterator VecR_it;
typedef VecR_vec::iterator VecR_it_non_const;

#endif
