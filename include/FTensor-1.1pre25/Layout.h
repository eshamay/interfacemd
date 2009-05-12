/* A simple enum to select what kind of storage order you want.
   Useful if you're mixing these classes with libraries with something
   that has row-major (not C/C++ standard) ordering like Fortran or
   some graphics libraries (Direct-X?). */

enum Layout {
  column_major,
  row_major
};
