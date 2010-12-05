#include "interface.h"

template <typename T>
double FindAveragedWaterInterface (T& t, const int numWaters, const bool top=true) {

	// load the waters and then sort them by position along the main axis
	t.LoadWaters();
	std::sort(t.int_wats.begin(), t.int_wats.end(), T::molecule_position_pred(Atom::O));

	// reverse the order if needed
	if (top)
		std::reverse (t.int_wats.begin(), t.int_wats.end());
	
	// grab the top numWaters of waters and find the average position of the along the axis
	VecR avg_pos (0.0,0.0,0.0);
	for (Mol_it it = t.int_wats.begin(); it != t.int_wats.begin() + numWaters; it++) {
		avg_pos += (*it)->GetAtom(Atom::O)->Position();
	}
	avg_pos = avg_pos * (1.0/numWaters);
	return avg_pos[t.axis];
}
