#ifndef GMXSYSTEM_H_
#define GMXSYSTEM_H_

#define CPLUSPLUS 1 

#include "mdsystem.h"
#include "trrfile.h"
#include "grofile.h"


typedef std::vector<Atom *> Atom_ptr_vec;
typedef std::vector<Molecule *> Mol_ptr_vec;

class GMXSystem : public MDSystem {

  public:
    GMXSystem (
	const char * trr_filepath,
	const char * gro_filepath);

    void LoadNext ();
    void LoadFirst ();

  private:
    TRRFile _trr;
    GROFile _gro;

    void _ParseMolecules ();
};

#endif
