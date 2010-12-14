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
#include <vector>
 
typedef enum {x=0, y=1, z=2} coord;

typedef Eigen::Vector3d VecR;
typedef Eigen::Vector3f VecF;
//typedef Eigen::MatrixBase<VecR>	vector_base;
typedef Eigen::MatrixBase<Eigen::Matrix<double, 3, 1, 2, 3, 1> >	vector_base;
//typedef Eigen::MapBase<Eigen::Matrix<double, 3, 1, 2, 3, 1> >		vector_map_base;
typedef Eigen::Map<VecR> vector_map;
typedef Eigen::MapBase<Eigen::Map<VecR, 1> >	vector_map_base;

typedef std::list< VecR, Eigen::aligned_allocator<VecR> > VecR_vec;
typedef VecR_vec VecR_list;
typedef VecR_vec::const_iterator VecR_it;
typedef VecR_vec::iterator VecR_it_non_const;

typedef std::list< VecF, Eigen::aligned_allocator<VecF> > VecF_vec;
typedef VecF_vec VecF_list;
typedef VecF_vec::const_iterator VecF_it;
typedef VecF_vec::iterator VecF_it_non_const;

// mapped coordinates
typedef Eigen::Map<VecR>	coord_t;
typedef std::vector<coord_t> coord_set_t;
typedef coord_set_t::const_iterator coord_it;

#endif
