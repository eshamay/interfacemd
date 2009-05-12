/* Takes a derivative of a Tensor0<T*> yielding a Tensor1.  It checks
   whether we're at a boundary, and takes a one-sided derivative
   there. */

template<class T, int Dim, char i>
class d_boundary_Tensor0
{
  const Tensor0<T*> &a;
  const Tensor1<int,Dim> &d_ijk;
  const Tensor1<double,Dim> &d_xyz;
  const Tensor2<bool,Dim,2> &boundary;
public:
  typename promote<T,double>::V operator()(const int N) const
  {
    return boundary(N,0) ? (*(&a+d_ijk(N))-a)*d_xyz(N) :
      (boundary(N,1) ? (a-*(&a-d_ijk(N)))*d_xyz(N) :
       (*(&a+d_ijk(N))-*(&a-d_ijk(N)))*d_xyz(N)*0.5);
  }
  d_boundary_Tensor0(const Tensor0<T*> &A,
		     const Tensor1<int,Dim> &D_ijk,
		     const Tensor1<double,Dim> &D_xyz,
		     const Tensor2<bool,Dim,2> &Boundary)
    : a(A), d_ijk(D_ijk), d_xyz(D_xyz), boundary(Boundary) {}
};

template<class T, int Dim, char i>
inline const Tensor1_Expr<const d_boundary_Tensor0<T,Dim,i>,typename promote<T,double>::V,Dim,i>
d_boundary(const Tensor0<T*> &a, const Index<i,Dim> index,
	   const Tensor1<int,Dim> &d_ijk,
	   const Tensor1<double,Dim> &d_xyz,
	   const Tensor2<bool,Dim,2> &boundary)
{
  typedef const d_boundary_Tensor0<T,Dim,i> Tensor_Expr;
  return Tensor1_Expr<Tensor_Expr,typename promote<T,double>::V,Dim,i>
    (Tensor_Expr(a,d_ijk,d_xyz,boundary));
}
