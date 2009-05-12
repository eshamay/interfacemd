/* A simple index class that tells the operators which template to
   instantiate.  The usual way to do this is to declare an index like

   Index<'i',3> i;

   It is important that each differently named index has a different
   template parameter. So

   Index<'i',3> i;
   Index<'j',3> j;          // Good

   is the right way, and

   Index<'i',3> i,j;        // Bad

   is likely to lead to errors, since the program thinks that i and j
   are identical. */

template<char i, int Dim>
class Index
{
public:
  Index() {};
};
