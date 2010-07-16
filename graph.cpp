#include "graph.h"

namespace bondgraph {
  BondGraph::Graph BondGraph::_graph (0);

  const double BondGraph::OHBONDLENGTH = 1.1;				// used to be 1.1
  const double BondGraph::HBONDLENGTH  = 2.3;				// used to be 2.46
  const double BondGraph::HBONDANGLECOS	= cos(30.0*M_PI/180.0);		// bonding angle has to be bigger than this cos (i.e. smaller than ~30 degrees
  const double BondGraph::NOBONDLENGTH = 2.0;
  const double BondGraph::NHBONDLENGTH = 1.3;		// uhmm... check this?
  const double BondGraph::SOBONDLENGTH = 1.8;



  // property map definitions
  BondGraph::PropertyMap<double,BondGraph::EdgeProperties>::Type BondGraph::b_length = get(&EdgeProperties::distance, _graph);
  BondGraph::PropertyMap<bondtype,BondGraph::EdgeProperties>::Type BondGraph::b_type = get(&EdgeProperties::btype, _graph);

  BondGraph::PropertyMap<Atom *,BondGraph::VertexProperties>::Type BondGraph::v_atom = get(&VertexProperties::atom, _graph);
  BondGraph::PropertyMap<VecR,BondGraph::VertexProperties>::Type BondGraph::v_position = get(&VertexProperties::position, _graph);
  BondGraph::PropertyMap<std::string,BondGraph::VertexProperties>::Type BondGraph::v_name = get(&VertexProperties::name, _graph);



  BondGraph::BondGraph ()
	:
	  _sys_type("xyz")
  { return; }



  BondGraph::BondGraph (const Atom_ptr_vec& atoms, const std::string& sys)
	:
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
	int i = 0;
	for (Atom_it it = atoms.begin(); it != atoms.end(); it++) {
	  v = boost::vertex(i, _graph);
	  v_atom[v] = *it;
	  v_position[v] = (*it)->Position();
	  v_name[v] = (*it)->Name();
	  i++;
	}

	return;
  }

  // some predicates
  bool BondGraph::_NameCombo (const std::string name1, const std::string name2, const std::string test1, const std::string test2) const {
	return 
	  ((name1.find(test1) != std::string::npos) || (name2.find(test1) != std::string::npos))
	  && 
	  ((name1.find(test2) != std::string::npos || name2.find(test2) != std::string::npos));
  }

  bool BondGraph::_SameAtomName (const std::string name1, const std::string name2) const {
	return (name1.find(name2) != std::string::npos) || (name2.find(name1) != std::string::npos);
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
		if (_SameAtomName(ai_name, aj_name)) continue;

		// calculate the distance between the two atoms
		double bondlength = MDSystem::Distance (v_position[*vi], v_position[*vj]).Magnitude();
		if (bondlength > HBONDLENGTH) continue;

		//double bondlength = v_position[*vi].MinDistance(v_position[*vj], Atom::_size);
		bondtype btype = unbonded;

		// first look at bonds between O and H
		if (_NameCombo(ai_name,aj_name,"O","H"))
		{
		  // one type of bond is the O-H covalent
		  if (bondlength <= OHBONDLENGTH) {
			btype = covalent;
		  }

		  // Or an H-bond is formed!
		  if (bondlength <= HBONDLENGTH && bondlength > OHBONDLENGTH) {
			// additionally, let's check the angle-criteria for an H-bond.
			// This is done by looking at the angle formed from
#ifdef ANGLE_CRITERIA
			// o1 is covalently bound to h, and o2 is h-bound to h
			Atom *o1 = (Atom *)NULL, *h = (Atom *)NULL, *o2 = (Atom *)NULL;	

			if (ai_name.find("O") != std::string::npos) {		// ai is the O, and aj is the H
			  o2 = ai;
			  h = aj;
			}
			else if (aj_name.find("O") != std::string::npos) {
			  h = ai;
			  o2 = aj;
			}

			o1 = h->ParentMolecule()->GetAtom("O");

			if (h == (Atom *)NULL || o1 == (Atom *)NULL || o2 == (Atom *)NULL) {
			  //throw (MALFORMED_H2O);
			  //printf ("Something wrong in assigning the atoms O and H in forming an H-bond - graph.cpp\n");
			  //exit(1);
			}

			VecR o1h = MDSystem::Distance (o1, h);	// the covalent bond
			VecR ho2 = MDSystem::Distance (h, o2);	// the H-bond

			double anglecos = o1h < ho2;		// cos(theta)

			if (anglecos > HBONDANGLECOS){
			  //printf ("% 10.3f / %6.3f\n", acos(angle)*180.0/M_PI, HBONDANGLECOS);
#endif
			  btype = hbond;
#ifdef ANGLE_CRITERIA
			}
#endif
		  }

		  /*
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
		   */
		}	// Check OH bond combos


		// now process SO2 molecules
		if (_NameCombo (ai_name,aj_name, "S", "O") && (bondlength < SOBONDLENGTH)) {
		  btype = covalent;
		}

		// add in the bond between two atoms
		this->_SetBond (*vi, *vj, bondlength, btype);
	  }
	}

	// Now fix up any weird atom-sharing between molecules. At this point we have to consider if we want to divide the system into separate molecules, or if we're interested in other phenomena, such as contact-ion pairs, etc.
	if (_sys_type == "xyz")
	  try {
		this->_ResolveSharedHydrogens ();
	  } catch (unboundhex& ex) {
		std::cout << "Exception caught while resolving hydrogens shared between multiple molecules" << std::endl;
		throw;
	  }
  }	// Parse Bonds

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
	for (Atom_it it = atoms.begin(); it != atoms.end(); it++)
	  boost::add_vertex(_graph);
	// parse the atom info into the vertices
	this->_ParseAtoms(atoms);
	// then find all the needed bond information
	try {
	  this->_ParseBonds ();
	} catch (graphex& ex) {
	  std::cout << "An exception was caught while determining the bonds between atoms in the system" << std::endl;
	}


	return;
  }

  void BondGraph::_SetBond (const Vertex& vi, const Vertex& vj, const double bondlength, const bondtype btype) {

	bool b;
	Edge e;

	tie(e, b) = add_edge(vi, vj, EdgeProperties(bondlength, btype), _graph);

	if (!b) {
	  /*
		 cout << "BondGraph::SetBond() - Tried to add a bond to an already-bonded atom pair" << endl;
		 v_atom[vi]->Print();
		 v_atom[vj]->Print();
		 exit(1);
	   */
	}

	return;
  }

  BondGraph::Vertex_it BondGraph::_FindVertex (Atom const * const ap) const {

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
		  if (v_name[*vi].find(name) == std::string::npos)
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
	for (Atom_it it = wat->begin(); it != wat->end(); it++) {
	  num += this->NumHBonds(*it);
	}

	return (num);
  }

  // calculates the water bonding coordination of a given water molecule
  coordination BondGraph::WaterCoordination (Water const * const wat) const
  {

	int c = 0;
	//Atom * ap;
	for (Atom_it atom = wat->begin(); atom != wat->end(); atom++) {
	  //for (int atom = 0; atom < wat->size(); atom++) {
	  //ap = wat->Atoms(atom);
	  //int bonds = NumHBonds (ap);
	  int bonds = NumHBonds (*atom);

	  if ((*atom)->Name().find("H") != std::string::npos) {
		c += 10 * bonds;
	  }
	  if ((*atom)->Name().find("O") != std::string::npos) {
		c += bonds;
	  }
	}

	return coordination(c);
  }	// water coordination


  void BondGraph::_ResolveSharedHydrogens () {

	// if any of the Hs are being shared between two molecules (i.e. a contact-ion pair) then we have to resolve which molecule gets the atom.
	// A simple solution is to decide that the molecule with the oxygen closer to the hydrogen is the one that wins out
	// This routine ***assumes*** that all the Hs will be bound to an Oxygen, and not some other atom in the system.
	// the closer oxygen becomes the covalent one, and the further oxygen becomes the h-bonded one

	// we'll go through each of the hydrogens and find if they are bound to one or two different molecules
	Vertex_it vi, vi_end, next;
	tie (vi, vi_end) = vertices(_graph);
	for (next = vi; vi != vi_end; vi = next) { next++;

	  Atom * H = _graph[*vi].atom;

	  if (H->Name().find("H") == std::string::npos) continue;

	  // here we check to see if it's bound to multiple atoms
	  Atom_ptr_vec covalent_atoms = this->BondedAtoms (H, covalent);

	  // great ... if the H is only covalently bound to one atom then don't worry about it
	  if (covalent_atoms.size() == 1) continue;

	  Atom_ptr_vec hbond_atoms = this->BondedAtoms (H, hbond);

	  // this is totally bizarre - since when do we see an H covalently bound to more than 1 molecule
	  //if (covalent_atoms.empty() > 2) {
	  /*
		 std::cout << "BondGraph::_ResolveSharedHydrogens()" << std::endl;
		 std::cout << "This H is covalently bound to more than 2 oxygens!!" << std::endl;
		 H->Print();
		 exit(1);
	   */
	  //throw (TOO_MANY_COVALENT_BONDS);
	  //}

	  // The hydrogen isn't close enough to any oxygen to be covalent - it's probably shared between 2 as hbonds
	  if (covalent_atoms.empty() && !hbond_atoms.empty()) {
		/*
		   std::cout << "\nBondGraph::_ResolveSharedHydrogens()" << std::endl;
		   std::cout << "Found an unbound H (no covalent bonds)!" << std::endl;
		   H->Print();
		   exit(1);
		 */

		// find the length of each of the H-bonds
		distance_vec distances;

		for (Atom_it it = hbond_atoms.begin(); it != hbond_atoms.end(); it++) {
		  distances.push_back (std::make_pair(this->Distance(H, *it), *it));
		}
		// sort nearest atoms
		md_utility::pair_sort_first (distances.begin(), distances.end());

		Vertex_it vO = this->_FindVertex(distances[0].second);

		// the closer bond is now changed from hbonded to covalent. The others are unchanged as they remain hbound
		this->_RemoveBond (*vi, *vO);								// cut the hbond to the nearest oxygen
		this->_SetBond(*vi, *vO, distances[0].first, covalent);	// and replace it with a covalent bond
	  }

	  // in the case that there are more than 1 covalent bond to the hydrogen (it's being shared, covalently, by more than two atoms!) then we have a physically very improbable scenario!
	  else if (covalent_atoms.size() >= 2) {
		std::cout << "Found a hydrogen atom that is covalently bound to multiple oxygens. This is physically very improbable." << std::endl;
		H->Print();
		throw (multiplyboundhex());
	  }

	  else if (covalent_atoms.empty() && hbond_atoms.empty()) {
		std::cout << "Found a hydrogen with no covalent or hydrogen bonds to which it is attached (below)" << std::endl;
		H->Print();
		throw (unboundhex());
	  }

	}

	return;
  }	// Resolve shared hydrogens

  double BondGraph::Distance (const Vertex& vi, const Vertex& vj) const {
	Edge e = edge(vi, vj, _graph).first;
	// problem!! if the second in the pair (bool) is false, then the edge doesn't exist! don't try to access edges between unbound atoms... or just add the edges to all the atom-pairs even if unbound. 
	if (!edge(vi, vj, _graph).second)  //----
	std::cout << b_length[e] << std::endl;

	return (b_length[e]);
  }

  double BondGraph::Distance (Atom const * const a1, Atom const * const a2) const {

	Vertex_it vi, vj;
	vi = this->_FindVertex (a1);
	vj = this->_FindVertex (a2);

	return Distance(*vi,*vj);
  }

  distance_pair BondGraph::ClosestAtom (const MolPtr& mol, const std::string& name) const {

	distance_vec distances;

	for (Atom_it it = mol->begin(); it != mol->end(); it++) {
	  distances.push_back (ClosestAtom (*it, name));
	}

	md_utility::pair_sort_first(distances.begin(), distances.end());

	return distances[0];

  }

  distance_pair BondGraph::ClosestAtom (const AtomPtr atom, const std::string& name) const {

	distance_vec distances;

	// find the distance between the target atom and all the other atoms in the system
	Vertex_it vt, vi, vi_end, next;
	vt = _FindVertex (atom);

	tie(vi, vi_end) = vertices(_graph);
	for (next = vi; vi != vi_end; vi = next) {
	  ++next;

	  // don't consider atoms within the same molecule, or an atom compared to itself
	  if (v_atom[*vt]->ParentMolecule() == v_atom[*vi]->ParentMolecule() || vt == vi) continue;

	  // do a name check so that only atoms with the (optional) given name are considered
	  if (name == "" || name == v_name[*vi])
		distances.push_back (std::make_pair (Distance(*vt, *vi), v_atom[*vi]));
	}

	// sort all the distances to find the one closest
	md_utility::pair_sort_first(distances.begin(), distances.end());

	return distances[0];
  }

  BondGraph::Edge BondGraph::_GetBond (const Vertex& vi, const Vertex& vj) const {

	Edge e;
	bool b;
	tie (e,b) = edge(vi, vj, _graph);

	if (!b) {
	  //throw(BOND_NOT_FOUND);
	  /*
		 printf ("BondGraph::_GetBond - Couldn't find the requested Bond");
		 exit(1);
	   */
	}

	return (e);
  }

  BondGraph::Edge BondGraph::_GetBond (Atom const * const a1, Atom const * const a2) const {

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

  }	// namespace bondgraph
