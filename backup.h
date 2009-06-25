void FindEulerAngles ();	// calculates the euler angles of the OH bonds treating each OH separately as the molec. Z-axis
	//double * EulerAngles (const int bond) { return _eulerAngles[bond]; }	// returns the euler angles orienting an OH-bond in the lab-frame


void Water::FindEulerAngles () {

/* method to determine the 3 euler angles in a space-fixed frame. This assumes that the space-fixed axes are described
   by the x, y, and z unit vectors given by [1,0,0], [0,1,0], and [0,0,1] respectively. The convention for the angle
   names and rotation directions is taken from "Angular Momentum: Understanding Spatial Aspects in Chemistry and Physics",
   Zare, Richard N; Wiley-Interscience 1988 (pp 77-81). This is sometimes referred to as the ZYZ ("y-") convention, with angles
   theta, phi, and chi.

   The procedure for determining the angles makes direct use of the direction cosine (rotation) matrix. The three angles are
   isolated by simple trigonometric and algebraic manipulation of the elements of the rotated vectors.

   The returned vector is composed of the three angles theta, phi, and chi, respectively.

   Note, we have to perform the procedure on both OH bonds simultaneously because of the molecule-rotation dependence of the psi angle
*/

	this->FindOHBonds();

	/* a quick thing that comes up with calculating spectra from slab systems is that the particles on one side are equivalent to the other side of the slab, but their spectra will cancel each other because of their opposite orientation... so let's flip the OH bonds for particles below the 1/2-way line of the system */
	/*
	if (_o->Z() < 26.0) {
		_oh1 *= -1.0;
		_oh2 *= -1.0;
	}
	*/

	// note: the molecular axes are defined where molAxes[0] is the molecular x-axis, molAxes[1] is the molecular y-axis, etc.
	// and molAxes[1][z] gives us the z-component of the molecular y-axis, and so on.

	// the first angle (theta) is the angle between the two Z axes. This is found from the direction cosine matrix quite directly
	_eulerAngles[0][0] = acos(_oh1.Unit()[z]);							// Theta for bond 1
	_eulerAngles[1][0] = acos(_oh2.Unit()[z]);							// Theta for bond 2

	// Phi is found similarly as a rotation about the lab-frame Z-axis. Also found from the cosine matrix.
	_eulerAngles[0][1] = atan2(_oh1.Unit()[y], _oh1.Unit()[x]);	// Phi of bond 1
	_eulerAngles[1][1] = atan2(_oh2.Unit()[y], _oh2.Unit()[x]);	// Phi of bond 2

	/* Now we have the theta and phi so we set up a rotation cosine matrix to transform from the lab-frame to the molecular frame.
	 * We still don't have psi, and can't have it quite yet, so we just rotate the theta and phi of both OH bonds, and then find psi from the rotated vectors
	 */

	// first we setup a rotation matrix from the lab from
	double rotate1[3][3];	// the rotation matrix
	double rotate2[3][3];	// the rotation matrix

	// we set up the angles initially only with the phi and theta we calculated, and leave psi as 0.0 - we calculate that later
	double angles1[3] = {_eulerAngles[0][0], _eulerAngles[0][1], 0.0};
	double angles2[3] = {_eulerAngles[1][0], _eulerAngles[1][1], 0.0};

	VecR rotatedOH1 = _oh1, rotatedOH2 = _oh2;	// the rotated vectors

	// We form the two rotation matrices based on the angle sets for both OH bonds
	// Then we rotate the OH vectors. Note the order here, we apply the rotation matrix 1 to OH vector 2, and vice versa
	//this->RotationFrom (rotatedOH2.Coords(), angles1);
	//this->RotationFrom (rotatedOH1.Coords(), angles2);

	// Now that we're in the molecular frame, we can calculate psi, the rotation about the molecular z-axis (like for phi above)
	_eulerAngles[0][2] = atan2(rotatedOH2[y], rotatedOH2[x]);	// Chi/Psi (depending on your choice of greek)
	_eulerAngles[1][2] = atan2(rotatedOH1[y], rotatedOH1[x]);

	// before leaving, let's restore the state of the OH bond back to what it was depending on location in the slab
	if (_o->Z() < 26.0) {
		_oh1 *= -1.0;
		_oh2 *= -1.0;
	}

return;
}
