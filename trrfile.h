#ifndef TRRFILE_H_
#define TRRFILE_H_

#define CPLUSPLUS 1

#include <iostream>
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_trr.h"
#include "vecr.h"

class TRRFile {

  public:
    TRRFile (char * path) : 
      _path(path),
      _file(xdrfile_open(path, "r"))
    { 
      read_trr_natoms(path, &_natoms);
      _x = new rvec[_natoms];
      _v = new rvec[_natoms];
      _f = new rvec[_natoms];

      _coords.resize(_natoms, VecR());
      _vels.resize(_natoms, VecR());
      _forces.resize(_natoms, VecR());

      this->LoadNext();

      return; 
    }

    ~TRRFile () {
      delete[] _x;
      delete[] _v;
      delete[] _f;
      xdrfile_close(_file);
      return;
    }

    int size () const { return _natoms; }

    void LoadFirst();
    void LoadNext();
    void PrintInfo () const;
    VecR_vec& Coords () { return _coords; }
    VecR& Coords (const int index) { return _coords[index]; }
    VecR& Velocities (const int index) { return _vels[index]; }
    VecR& Forces (const int index) { return _forces[index]; }

  private:

    char * _path;
    XDRFILE * _file;
    int	_natoms;
    int _frame;
    float _time;
    float _lambda;
    matrix _box;
    rvec *_x, *_v, *_f;
    VecR_vec _coords, _vels, _forces;

    void PrintBox () const;
    void Print (const VecR_vec& vec) const;
};
#endif
