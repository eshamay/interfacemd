#pragma once
#ifndef GROFILE_H_
#define GROFILE_H_

#include "atom.h"
#include "molecule.h"
#include <string>
#include <vector>
#include <iostream>

class GROFile {
  public:
    GROFile (const std::string path) :
      _file (fopen (path.c_str(), "r"))
    { 
      if (_file == (FILE *)NULL) {
	std::cout << "Error opening the .gro file " << path << std::endl;
	exit(1);
      }
      _ParseSystem();
      return; 
    }

    ~GROFile () 
    {
      fclose(_file);
      for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++)
	delete ((*it));
      for (Mol_it it = _mols.begin(); it != _mols.end(); it++)
	delete ((*it));

      return;
    }

    void Print () const;
    int NumAtoms () const { return _natoms; }
    Atom_ptr_vec& Atoms() { return _atoms; }
    Mol_ptr_vec& Molecules() { return _mols; }

    AtomPtr operator[] (const int index) { return _atoms[index]; }

  private:
    FILE * _file;
    std::string _title;
    int _natoms;
    VecR _box_size;

    Atom_ptr_vec _atoms;
    Mol_ptr_vec _mols;

    void _ParseSystem ();
};

#endif
