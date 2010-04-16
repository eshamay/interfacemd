#include "trrfile.h"

void TRRFile::LoadFirst () {
  xdrfile_close(_file);
  _file = xdrfile_open(_path, "r");
  this->LoadNext();

  return;
}

// load the next frame of the trr file - coordinates, velocities, forces
void TRRFile::LoadNext () {
  // load the frame using the xdrfile framework
  read_trr(
      _file, _natoms, &_frame, &_time, &_lambda, _box, _x, _v, _f);

  for (int i = 0; i < _natoms; i++) {
    _coords[i].Set(_x[i][0], _x[i][1], _x[i][2]);
    _vels[i].Set(_v[i][0], _v[i][1], _v[i][2]);
    _forces[i].Set(_f[i][0], _f[i][1], _f[i][2]);
  }

 return;
}

void TRRFile::PrintInfo () const {
  printf ("Number of atoms = %d\n", _natoms);
  printf ("Frame #%d\n", _frame);
  printf ("Time = % 8.3f ps\n", _time);
  printf ("Lambda = % 8.3f\n", _lambda);

  printf ("\nBox:\n");
  this->PrintBox();

  printf ("\nCoordinates: \n");
  this->Print(_coords);
  printf ("\n");
 
  printf ("\nVelocities: \n");
  this->Print(_vels);
  printf ("\n");

  printf ("\nForces: \n");
  this->Print(_forces);
  printf ("\n");

  return;
}

void TRRFile::PrintBox () const {

  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      printf ("% 13.4f", _box[i][j]);
    }
    printf ("\n");
  }

  return;
}

void TRRFile::Print (const VecR_vec& vec) const {
  for (VecR_it it = vec.begin(); it != vec.end(); it++) {
    it->Print();
  }
  printf ("\n");
  return;
}

/*
int main () {
  TRRFile trr ("sw_md.trr");
  trr.PrintInfo();
  trr.LoadNext();
  trr.PrintInfo();
  return 0;
}
*/
