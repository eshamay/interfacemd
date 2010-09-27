#pragma once
#ifndef VECR_H_
#define VECR_H_

#ifndef EIGEN_MATRIXBASE_PLUGIN 
#define EIGEN_MATRIXBASE_PLUGIN "EigenMatrixAddon.h"
#endif

#include <Eigen/Core>
#include <Eigen/Geometry>
USING_PART_OF_NAMESPACE_EIGEN

#include <list>
 
typedef enum {x=0, y=1, z=2} coord;

typedef Eigen::Vector3d VecR;
typedef Eigen::Vector3f VecF;

typedef std::list< VecR, Eigen::aligned_allocator<VecR> > VecR_vec;
typedef VecR_vec VecR_list;
typedef VecR_vec::const_iterator VecR_it;
typedef VecR_vec::iterator VecR_it_non_const;

typedef std::list< VecF, Eigen::aligned_allocator<VecF> > VecF_vec;
typedef VecF_vec VecF_list;
typedef VecF_vec::const_iterator VecF_it;
typedef VecF_vec::iterator VecF_it_non_const;

#endif
