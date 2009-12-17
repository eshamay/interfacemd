#include "graph.h"

Graph BondGraph::_graph (0);

// property map definitions
PropertyMap<double,EdgeProperties>::Type BondGraph::b_length = get(&EdgeProperties::distance, _graph);
PropertyMap<bondtype,EdgeProperties>::Type BondGraph::b_type = get(&EdgeProperties::btype, _graph);

PropertyMap<Atom *,VertexProperties>::Type BondGraph::v_atom = get(&VertexProperties::atom, _graph);
PropertyMap<VecR,VertexProperties>::Type BondGraph::v_position = get(&VertexProperties::position, _graph);
PropertyMap<std::string,VertexProperties>::Type BondGraph::v_name = get(&VertexProperties::name, _graph);

BondGraph::BondGraph ()
{ return; }

BondGraph::BondGraph (const Atom_ptr_vec& atoms, std::string sys) :
	_sys_type(sys)
{

	this->UpdateGraph(atoms);

return;
}

BondGraph::~BondGraph () {

	this->_ClearBonds();
	this->_ClearAtoms();

return;
}

void BondGraph::_ParseAtoms (const Atom_ptr_vec& atoms) {

	// set all the vertices to contain the proper info
	Vertex v;
	RUN (atoms) {
		v = boost::vertex(i, _graph);
		v_atom[v] = atoms[i];
		v_position[v] = atoms[i]->Position();
		v_name[v] = atoms[i]->Name();
	}

return;
}

// Calculate the distance between each of the atoms and update the graph edge list with bonds between them
// Edges/bonds should only exist between atoms that are bound in some fashion (h-hond, covalent, etc.)
void BondGraph::_ParseBonds () {

  // first clear out all the bonds from before
  this->_ClearBonds();

#ifdef ANGLE_CRITERIA
  // now run through all atom combos and get their bondlengths in order to set their bondtypes
  Atom *ai, *aj;
#endif

  // this little for-loop bit runs through all atom-pair combinations once and find the bond-types between them
  Vertex_it vi, vi_end, next_i;
  tie(vi, vi_end) = vertices(_graph);
  Vertex_it vj, vj_end, next_j;
  tie(vj, vj_end) = vertices(_graph);

  --vi_end;
  for (next_i = vi; vi != vi_end; vi = next_i) {
	++next_i;
	vj = next_i;
	for (next_j = vj; vj != vj_end; vj = next_j) {
	  ++next_j;

#ifdef ANGLE_CRITERIA
	  ai = v_atom[*vi];
	  aj = v_atom[*vj];
#endif
	  std::string ai_name = v_name[*vi];
	  std::string aj_name = v_name[*vj];

	  // Don't connect oxygens to oxygens, and hydrogen to hydrogen...etc.
	  if (ai_name == aj_name) continue;

	  // calculate the distance between the two atoms
	  double bondlength = MDSystem::Distance (v_position[*vi], v_position[*vj]).Magnitude();
	  //double bondlength = v_position[*vi].MinDistance(v_position[*vj], Atom::_size);
	  bondtype btype = unbonded;

	  // first look at bonds between O and H
	  if ( (ai_name.find("O") != std::string::npos || aj_name.find("O") != std::string::npos)
		  &&
		  (ai_name.find("H") != std::string::npos || aj_name.find("H") != std::string::npos)
		 )
	  {
		// one type of bond is the O-H covalent
		if (bondlength <= OHBONDLENGTH) {
		  btype = ohbond;
		}

		// Or an H-bond is formed!
		if (bondlength <= HBONDLENGTH && bondlength > OHBONDLENGTH) {
		  // additionally, let's check the angle-criteria for an H-bond.
		  // This is done by looking at the angle formed from
#ifdef ANGLE_CRITERIA
		  Atom *o1, *h, *o2;	// o1 is covalently bound to h, and o2 is h-bound to h
		  Water * wat;

		  if (ai_name.find("O") != std::string::npos) {		// ai is the O, and aj is the H
			o2 = ai;
			h = aj;
		  }
		  else if (aj_name.find("O") != std::string::npos) {
			h = ai;
			o2 = aj;
		  }
		  wat = static_cast<Water *>(h->ParentMolecule());
		  o1 = wat->GetAtom("O");

		  VecR oh1 = MDSystem::Distance (o1, h);	// the covalent bond
		  VecR oh2 = MDSystem::Distance (o2, h);	// the H-bond

		  double angle = acos(oh1 < oh2);
		  //printf ("% 10.3f\n", angle * 180.0/M_PI);

		  if (angle < HBONDANGLE)
#endif

			btype = hbond;
		}
	  }

	  // now connect Os to Ns
	  if ( (ai_name.find("O") != std::string::npos || aj_name.find("O") != std::string::npos)
		  &&
		  (ai_name.find("N") != std::string::npos || aj_name.find("N") != std::string::npos)
		 )
	  {
		if (bondlength <= NOBONDLENGTH) {
		  btype = nobond;
		}
	  }
	  // add in the bond between two atoms
	  if (btype != unbonded)
		this->_SetBond (*vi, *vj, bondlength, btype);
	}}

	// Now fix up any weird atom-sharing between molecules. At this point we have to consider if we want to divide the system into separate molecules, or if we're interested in other phenomena, such as contact-ion pairs, etc.
	if (_sys_type == "xyz")
	  this->_FixSharedAtoms ();
	return;
}

// Set all the bonds to unbonded
void BondGraph::_ClearBonds () {
  // Remove all the edges.
  Edge_it ei, e_end, next;
  tie(ei, e_end) = edges(_graph);
  for (next = ei; ei != e_end; ei = next) {
	++next;
		remove_edge(*ei, _graph);
	}
return;
}

void BondGraph::_ClearAtoms () {
	// Remove all the vertices.
	Vertex_it vi, vi_end, next;
	tie(vi, vi_end) = vertices(_graph);
	for (next = vi; vi != vi_end; vi = next) {
		++next;
		remove_vertex(*vi, _graph);
	}
return;
}

void BondGraph::UpdateGraph (const Atom_ptr_vec& atoms) {

	// clear out the old graph info
	_graph.clear();
	// resize the graph with new vertices
	RUN (atoms)
		boost::add_vertex(_graph);
	// parse the atom info into the vertices
	this->_ParseAtoms(atoms);
	// then find all the needed bond information
	this->_ParseBonds ();


return;
}

void BondGraph::_SetBond (const Vertex& vi, const Vertex& vj, const double bondlength, const bondtype btype) {

	bool b;
	Edge e;

	tie(e, b) = add_edge(vi, vj, EdgeProperties(bondlength, btype), _graph);

	if (!b) {
		cout << "BondGraph::SetBond() - Tried to add a bond to an already-bonded atom pair" << endl;
		v_atom[vi]->Print();
		v_atom[vj]->Print();
		exit(1);
	}

return;
}

Vertex_it BondGraph::_FindVertex (Atom const * const ap) const {

	Vertex_it vi, vi_end, next;
	tie(vi, vi_end) = vertices(_graph);
	for (next = vi; vi != vi_end; vi = next) {
		++next;
		if (v_atom[*vi] == ap)
			break;
	}

return (vi);
}

// returns a list of the atoms bonded to the given atom
Atom_ptr_vec BondGraph::BondedAtoms
	(Atom const * const ap,
	bondtype const btype,
	std::string const name)
const {

	Atom_ptr_vec atoms;
	Vertex_it va = this->_FindVertex(ap);

	Adj_it vi, vi_end, next;
	tie(vi, vi_end) = adjacent_vertices(*va, _graph);
	for (next = vi; vi != vi_end; vi = next) {
		next++;
		// check the bondtype criteria
		if (btype != unbonded) {
			Edge e = edge(*va, *vi, _graph).first;

			if (
				(btype == covalent && b_type[e] == hbond)
				||
				(btype != covalent && b_type[e] != btype))
					continue;
		}

		// check the name criteria
		if (name != "") {
			if (v_name[*vi].find(name) == string::npos)
				continue;
		}

		atoms.push_back(v_atom[*vi]);
	}

return (atoms);
}

int BondGraph::NumHBonds (Atom const * const ap) const {

	return (this->BondedAtoms(ap, hbond).size());
}

int BondGraph::NumHBonds (Water const * const wat) const {

	int num = 0;
	RUN (wat->Atoms()) {
		num += this->NumHBonds(wat->Atoms(i));
	}

return (num);
}

// calculates the water bonding coordination of a given water molecule
coordination BondGraph::WaterCoordination (Water const * const wat) const {

	int c = 0;
	Atom * ap;
	for (int atom = 0; atom < wat->size(); atom++) {
		ap = wat->Atoms(atom);
		int bonds = NumHBonds (ap);

		if (ap->Name().find("H") != std::string::npos) {
			c += 10 * bonds;
		}
		if (ap->Name().find("O") != std::string::npos) {
			c += bonds;
		}
	}
//	printf ("%d)  ", c);

return (coordination) c;
}

void BondGraph::_FixSharedAtoms () {

		// if any of the Hs are being shared between two molecules (i.e. a contact-ion pair) then we have to resolve which molecule gets the atom.
		// A simple solution is to decide that the molecule with the oxygen closer to the hydrogen is the one that wins out
		// This routine ***assumes*** that all the Hs will be bound to an Oxygen, and not some other atom in the system.

	// we'll go through each of the hydrogens and find if they are bound to one or two different molecules
	Vertex_it vi, vi_end, next;
	tie (vi, vi_end) = vertices(_graph);
	for (next = vi; vi != vi_end; vi = next) { next++;

		Atom * H = _graph[*vi].atom;

		if (H->Name().find("H") == string::npos) continue;

		// here we check to see if it's bound to multiple atoms
		std::vector<Atom *> atoms = this->BondedAtoms (H, covalent);

		// great ... if the H is only bound to one molecule
		if (atoms.size() == 1) continue;

		// we shouldn't have any free-floating hydrogens!
		if (atoms.size() < 1) {
			std::cout << "BondGraph::_FixSharedAtoms()" << std::endl;
			std::cout << "Found an unbound H!" << std::endl;
			H->Print();
			exit(1);
		}

		// this is totally bizarre - since when do we see an H bound to 3 molecules
		if (atoms.size() > 2) {
			std::cout << "BondGraph::_FixSharedAtoms()" << std::endl;
			std::cout << "This H is bound to more than 2 molecules!!" << std::endl;
			H->Print();
			exit(1);
		}

		// now we do a comparison of distances between the two covalently bound atoms and the hydrogen - which is further away?
		int further;
		double distance[2] = { this->Distance(H, atoms[0]), this->Distance(H, atoms[1]) };

		further = (distance[0] < distance[1]) ? 1 : 0;

		Vertex_it vO = this->_FindVertex(atoms[further]);

		// the further bond is now changed from covalent to hbonded just to show that it's not the primary bond to the hydrogen
		this->_RemoveBond (*vi, *vO);
		this->_SetBond(*vi, *vO, distance[further], hbond);
	}

return;
}

double BondGraph::Distance (Atom const * const a1, Atom const * const a2) const {

	Vertex_it vi, vj;
	vi = this->_FindVertex (a1);
	vj = this->_FindVertex (a2);

	Edge e = edge(*vi, *vj, _graph).first;

return (b_length[e]);
}

Edge BondGraph::_GetBond (const Vertex& vi, const Vertex& vj) const {

	Edge e;
	bool b;
	tie (e,b) = edge(vi, vj, _graph);

	if (!b) {
		printf ("BondGraph::_GetBond - Couldn't find the requested Bond");
		exit(1);
	}

return (e);
}

Edge BondGraph::_GetBond (Atom const * const a1, Atom const * const a2) const {

	Vertex_it v1 = this->_FindVertex(a1);
	Vertex_it v2 = this->_FindVertex(a2);
	Edge e = this->_GetBond(*v1, *v2);

return (e);
}

void BondGraph::_RemoveBond (const Vertex& vi, const Vertex& vj) {

	Edge e = this->_GetBond(vi, vj);
	remove_edge(e, _graph);

return;
}

void BondGraph::_RemoveBond (Atom const * const a1, Atom const * const a2) {

	Edge e = this->_GetBond (a1, a2);
	remove_edge(e, _graph);

return;
}
