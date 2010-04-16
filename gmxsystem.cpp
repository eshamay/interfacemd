#include "gmxsystem.h"

GMXSystem::GMXSystem (const char * trr_filepath, const char * gro_filepath)
:
  _trr(const_cast<char *>(trr_filepath)), _gro(gro_filepath)
{ 
  this->LoadFirst();
  return; 
}

void GMXSystem::LoadFirst() {
  _trr.LoadFirst();
  _mols = _gro.Molecules();
  _atoms = _gro.Atoms();

  this->_ParseMolecules();

  return;
}

void GMXSystem::LoadNext() {
  _trr.LoadNext();
  this->_ParseMolecules();
}

void GMXSystem::_ParseMolecules () {
  for (int i = 0; i < _gro.NumAtoms(); i++) {
    _atoms[i]->Position(_trr.Coords(i));
    _atoms[i]->Force(_trr.Forces(i));
  }
}

/*
int main () {

  GMXSystem sys ("sw_md.trr", "sw_md.gro");
  Molecule * mol = sys.Molecules(0);
  for (int i = 0; i < 10; i++) {
    mol->Print();
    sys.LoadNext();
  }

  return 0;
}
*/
