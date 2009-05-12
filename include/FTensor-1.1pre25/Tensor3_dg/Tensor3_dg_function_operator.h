/* This contains the definitions of the almost all of the indexing
   operators for Tensor3_dg. */

/* These operator()'s are the first part in constructing template
   expressions. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, char k, int Dim01, int Dim2>
inline Tensor3_dg_Expr<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
  T,Dim01,Dim2,i,j,k>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()
  (const Index<i,Dim01> index1, const Index<j,Dim01> index2,
   const Index<k,Dim2> index3)
{
  return Tensor3_dg_Expr<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
    T,Dim01,Dim2,i,j,k>(*this);
}

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, char k, int Dim01, int Dim2>
inline Tensor3_dg_Expr<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
  T,Dim01,Dim2,i,j,k>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()
  (const Index<i,Dim01> index1, const Index<j,Dim01> index2,
   const Index<k,Dim2> index3) const
{
  return Tensor3_dg_Expr<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
    T,Dim01,Dim2,i,j,k> (*this);
}

/* These operators are for internal contractions. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, int Dim, int Dim12>
inline Tensor1_Expr<const Tensor3_contracted_12
<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,Dim12,i>,T,Dim,i>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()
  (const Index<i,Dim> index1, const Index<j,Dim12> index2,
   const Index<j,Dim12> index3) const
{
  typedef const Tensor3_contracted_12<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
    T,Dim12,i> TensorExpr;
  return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this));
}

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, int Dim, int Dim02>
inline Tensor1_Expr<const Tensor3_contracted_02
<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,Dim02,i>,T,Dim,i>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()
  (const Index<j,Dim02> index1, const Index<i,Dim> index2,
   const Index<j,Dim02> index3) const
{
  typedef const Tensor3_contracted_02<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
    T,Dim02,i> TensorExpr;
  return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this));
}

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, int Dim, int Dim01>
inline Tensor1_Expr<const Tensor3_contracted_01
<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,Dim01,i>,T,Dim,i>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()
  (const Index<j,Dim01> index1, const Index<j,Dim01> index2,
   const Index<i,Dim> index3) const
{
  typedef const Tensor3_contracted_01<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
    T,Dim01,i> TensorExpr;
  return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this));
}

/* This is for expressions where a number is used for one slot, and
   indices for the others, yielding a Tensor2_Expr or
   Tensor2_symmetric_Expr.  The non-const versions don't actually
   create a Tensor3_dg_number_rhs_* object, while the const versions
   do create a Tensor3_dg_number_*. */

/* First slot. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, int N, int Dim0, int Dim1>
Tensor2_Expr<Tensor3_dg_number_rhs_0<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N>,
  T,Dim0,Dim1,i,j>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Number<N> n1,
				     const Index<i,Dim0> index1,
				     const Index<j,Dim1> index2)
{
  typedef Tensor3_dg_number_rhs_0<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N> TensorExpr;
  return Tensor2_Expr<TensorExpr,T,Dim0,Dim1,i,j>(*this);
}

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, int N, int Dim0, int Dim1>
const Tensor2_Expr<const Tensor3_dg_number_0<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N>,
  T,Dim0,Dim1,i,j>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Number<N> n1,
				     const Index<i,Dim0> index1,
				     const Index<j,Dim1> index2) const
{
  typedef const Tensor3_dg_number_0<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N> TensorExpr;
  return Tensor2_Expr<TensorExpr,T,Dim0,Dim1,i,j>(TensorExpr(*this));
}

/* Second slot. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, int N, int Dim0, int Dim1>
Tensor2_Expr<Tensor3_dg_number_rhs_0<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N>,
  T,Dim0,Dim1,i,j>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Index<i,Dim0> index1,
				     const Number<N> n1,
				     const Index<j,Dim1> index2)
{
  typedef Tensor3_dg_number_rhs_0<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N> TensorExpr;
  return Tensor2_Expr<TensorExpr,T,Dim0,Dim1,i,j>(*this);
}

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, int N, int Dim0, int Dim1>
const Tensor2_Expr<const Tensor3_dg_number_0<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N>,
  T,Dim0,Dim1,i,j>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Index<i,Dim0> index1,
				     const Number<N> n1,
				     const Index<j,Dim1> index2) const
{
  typedef const Tensor3_dg_number_0<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N> TensorExpr;
  return Tensor2_Expr<TensorExpr,T,Dim0,Dim1,i,j>(TensorExpr(*this));
}

/* Third slot. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, int N, int Dim>
inline Tensor2_symmetric_Expr<Tensor3_dg_number_rhs_2<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N>,
  T,Dim,i,j>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Index<i,Dim> index1,
				     const Index<j,Dim> index2,
				     const Number<N> n1)
{
  typedef Tensor3_dg_number_rhs_2<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,T,Dim,i,j>(*this);
}

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, int N, int Dim>
const Tensor2_symmetric_Expr<const Tensor3_dg_number_2
<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N>,T,Dim,i,j>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Index<i,Dim> index1,
				     const Index<j,Dim> index2,
				     const Number<N> n1) const
{
  typedef const Tensor3_dg_number_2<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,T,Dim,i,j>(TensorExpr(*this));
}

/* This is for expressions where a number is used for two slots, and
   an Index for the other, yielding a Tensor1_Expr.  The non-const
   versions don't actually create a Tensor3_dg_number_rhs_* object,
   while the const versions do create a Tensor3_dg_number_*. */

/* Index in first slot. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, int N1, int N2, int Dim>
Tensor1_Expr<Tensor3_dg_number_rhs_12<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N1,N2>,T,Dim,i>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Index<i,Dim> index,
				     const Number<N1> n1,
				     const Number<N2> n2)
{
  typedef Tensor3_dg_number_rhs_12<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N1,N2> TensorExpr;
  return Tensor1_Expr<TensorExpr,T,Dim,i>(*this);
}

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, int N1, int N2, int Dim>
const Tensor1_Expr<const Tensor3_dg_number_12<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
  T,N1,N2>,T,Dim,i>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Index<i,Dim> index,
				     const Number<N1> n1,
				     const Number<N2> n2) const
{
  typedef const Tensor3_dg_number_12<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N1,N2>
    TensorExpr;
  return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this));
}

/* Index in second slot.  I use the same structures as for the Index
   in the first slot since the tensor is symmetric on the first two
   indices. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, int N1, int N2, int Dim>
Tensor1_Expr<Tensor3_dg_number_rhs_12<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N1,N2>,T,Dim,i>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Number<N1> n1,
				     const Index<i,Dim> index,
				     const Number<N2> n2)
{
  typedef Tensor3_dg_number_rhs_12<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N1,N2> TensorExpr;
  return Tensor1_Expr<TensorExpr,T,Dim,i>(*this);
}

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, int N1, int N2, int Dim>
const Tensor1_Expr<const Tensor3_dg_number_12<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
  T,N1,N2>,T,Dim,i>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Number<N1> n1,
				     const Index<i,Dim> index,
				     const Number<N2> n2) const
{
  typedef const Tensor3_dg_number_12<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N1,N2>
    TensorExpr;
  return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this));
}

/* Index in third slot. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, int N1, int N2, int Dim>
Tensor1_Expr<Tensor3_dg_number_rhs_01<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N1,N2>,T,Dim,i>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Number<N1> n1, const Number<N2> n2,
				     const Index<i,Dim> index)
{
  typedef Tensor3_dg_number_rhs_01<Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N1,N2> TensorExpr;
  return Tensor1_Expr<TensorExpr,T,Dim,i>(*this);
}

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, int N1, int N2, int Dim>
const Tensor1_Expr<const Tensor3_dg_number_01<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
  T,N1,N2>,T,Dim,i>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Number<N1> n1, const Number<N2> n2,
				     const Index<i,Dim> index) const
{
  typedef const Tensor3_dg_number_01<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T,N1,N2>
    TensorExpr;
  return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this));
}

/* Specializations for using actual numbers instead of Number<>.
   This is for expressions where an actual number is used for one
   slot, and indices for the others, yielding a Tensor2_Expr or
   Tensor2_symmetric_Expr. I only define the const versions. */

/* First slot. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, int Dim0, int Dim1>
const Tensor2_Expr<const Tensor3_dg_numeral_0<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T>,
  T,Dim0,Dim1,i,j>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const int N, const Index<i,Dim0> index1,
				     const Index<j,Dim1> index2) const
{
  typedef const Tensor3_dg_numeral_0<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T> TensorExpr;
  return Tensor2_Expr<TensorExpr,T,Dim0,Dim1,i,j>(TensorExpr(*this,N));
}

/* Second slot. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, int Dim0, int Dim1>
const Tensor2_Expr<const Tensor3_dg_numeral_0<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T>,
  T,Dim0,Dim1,i,j>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Index<i,Dim0> index1, const int N,
				     const Index<j,Dim1> index2) const
{
  typedef const Tensor3_dg_numeral_0<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T> TensorExpr;
  return Tensor2_Expr<TensorExpr,T,Dim0,Dim1,i,j>(TensorExpr(*this,N));
}

/* Third slot. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, char j, int Dim>
const Tensor2_symmetric_Expr<const Tensor3_dg_numeral_2
<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T>,T,Dim,i,j>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Index<i,Dim> index1,
				     const Index<j,Dim> index2,
				     const int N) const
{
  typedef const Tensor3_dg_numeral_2<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,T,Dim,i,j>(TensorExpr(*this,N));
}

/* This is for expressions where an actual number is used for two
   slots, and an Index for the other, yielding a Tensor1_Expr. I
   only define the const versions. */

/* Index in first slot. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, int Dim>
const Tensor1_Expr<const Tensor3_dg_numeral_12<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
  T>,T,Dim,i>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const Index<i,Dim> index, const int N1,
				     const int N2) const
{
  typedef const Tensor3_dg_numeral_12<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T>
    TensorExpr;
  return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this,N1,N2));
}

/* Index in second slot.  I use the same structures as for the Index
   in the first slot since the tensor is symmetric on the first two
   indices. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, int Dim>
const Tensor1_Expr<const Tensor3_dg_numeral_12<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
  T>,T,Dim,i>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const int N1, const Index<i,Dim> index,
				     const int N2) const
{
  typedef const Tensor3_dg_numeral_12<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T>
    TensorExpr;
  return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this,N1,N2));
}

/* Index in third slot. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
template<char i, int Dim>
const Tensor1_Expr<const Tensor3_dg_numeral_01<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,
  T>,T,Dim,i>
Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>::operator()(const int N1, const int N2,
				     const Index<i,Dim> index) const
{
  typedef const Tensor3_dg_numeral_01<const Tensor3_dg<T,Tensor_Dim01,Tensor_Dim2>,T>
    TensorExpr;
  return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this,N1,N2));
}
