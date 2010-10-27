#include "watersfg.h"

namespace morita {
	double SFGCalculator::alpha_data1[9] = {1.539, 0.000, -0.163, 0.000, 1.656, 0.000, -0.163, 0.000, 7.200};
	double SFGCalculator::alpha_data2[9] = {6.4684, 0.0, -1.9057, 0.0, 1.656, 0.0, -1.9057, 0.0, 2.2703};	// this is from a unitary transformation, oh-oh angle of 109.4719

	// this is for a direct rotation D*alpha, without any unitary transform business taking place. Just direct left-multiplication by the DCM from oh2 to oh1 using the POL3 109.4719 angle.
	//double alpha_data2[9] = {0.3588, 0.0, 6.7339, 0.0, -1.656, 0.0, 1.5053, 0.0, -2.5534};
	MatR SFGCalculator::AlphaDeriv1 (alpha_data1);
	MatR SFGCalculator::AlphaDeriv2 (alpha_data2);

	VecR SFGCalculator::MuDeriv1 (-0.058, 0.000, 0.157);	// dipole derivative vector from the paper
	VecR SFGCalculator::MuDeriv2 (0.1287, 0.0, -0.1070);	// dipole derivative of the 2nd OH bond, in the frame of the first (found by direction cosine rotation)

	//SFGCalculator::SFGCalculator (string const polarization, coord const axis) {
	SFGCalculator::SFGCalculator () { }

	/* taken from Eq. 10c from the Morita-Hynes work to calculate the change in frequency of a bond under a force. Returns two frequency shifts due to forces on a molecule */
	void SFGCalculator::FreqShift (morita::MoritaH2O& water) {

		Atom * o_atom = water.O();
		Atom * h1_atom = water.H1();
		Atom * h2_atom = water.H2();

		// first we grab the force vectors on each of the atoms	(Units = kcal/mol/Angstrom)
		VecR vForceO (o_atom->Force());
		VecR vForceH1 (h1_atom->Force());
		VecR vForceH2 (h2_atom->Force());

		/******************************
		 * Calculate the forces on the two bonds
		 ******************************/

		/* The force calculation - splitting up the forces on the atoms between the two bonds - is done by find the inner-product of the force with the OH bond. Another method to try is also to scale the oxygen force and hydrogen force based on their distance from center of mass, angle ratios of the force vector with the oh bond, or based on mass of the atom as a percentage of the total mass, etc.
		*/

		// we need certain values before being able to proceed. We'll need the two target water OH's. And in order to generate those, we'll need the target water geometry. Let's start with the OH bond vectors (pointing from O to H).
		VecR const * const oh1 = water.OH1();
		VecR const * const oh2 = water.OH2();
		// and the location of the target water's oxygen
		VecR o (o_atom->Position());

		// the forces on the hydrogens are scaled by the distance from the center of mass to the hydrogen and are taken as the dot product of the force vector with that scaled bond vector.
		//double Hscale = 0.6;
		double Hscale = (OH_LENGTH-OH_COM_LENGTH)/sfg_units::ANG2BOHR;
		double ForceH1 = vForceH1.dot((oh1->normalized() * Hscale));
		double ForceH2 = vForceH2.dot((oh2->normalized() * Hscale));

		// here are the angles between the oxygen force vectors and the vectors pointing from the oxygens to the center of masses
		// these below are the cosines of the angles
		double cos1 = (vForceO < oh1->normalized());
		double cos2 = (vForceO < oh2->normalized());

		// to get the oxygen force distributed we do a few things:
		// 	first we need to determine the amount of the oxygen force that goes into each bond. One way is to find some ratio of angles between the force vector and the bonds. if we have theta1 as the angle between the force and OH bond1, and theta2 for the other, then we should be able to work a way to distribute the force completely.
		// 		cos1/fabs(cos1) gives us a sense of the direction of the force - is it compressing or stretching the bond. Negative values imply that the bond is being stretched. Because we're interested in the "stretching" force on the bond, we take the negative of this ratio.
		// 		The angle ratio is tempramental... still working that out
		double Oscale = OH_COM_LENGTH/sfg_units::ANG2BOHR;
		//double Oscale = 0.8;
		double ForceO1 = (vForceO.dot((oh1->normalized() * Oscale))) * fabs(cos1)/(fabs(cos1)+fabs(cos2));
		double ForceO2 = (vForceO.dot((oh2->normalized() * Oscale))) * fabs(cos2)/(fabs(cos1)+fabs(cos2));

		double ForceOH1 = ForceH1 - ForceO1;
		double ForceOH2 = ForceH2 - ForceO2;
		//printf ("% 10.3f\t% 10.3f\n", ForceOH1, ForceOH2);

		// convert the forces into atomic units
		ForceOH1 *= sfg_units::AMBER2ATOMIC;
		ForceOH2 *= sfg_units::AMBER2ATOMIC;

		/******************************
		 * Calculate frequency shifts due to forces compressing or stretching bonds
		 ******************************/
		/* now perform the direct application of the equation Eq.10c */
		// note: the frequency shift is given as an angular frequency (omega) in the equation. Do we need a conversion to Hz style frequency here?
		_OH1FreqShift = PREFACTOR * ForceOH1;		// now in frequency in atomic units
		_OH2FreqShift = PREFACTOR * ForceOH2;

		//printf ("% 10.3f\t% 10.3f\n", _OH1FreqShift*AU2WAVENUMBER, _OH2FreqShift*AU2WAVENUMBER);

#ifdef	DIPOLE_DIPOLE		// only if adding in the dipole-dipole correction term (which should be pretty small)
		std::vector<double> dipolePotential (2, 0.0);
		this->DipoleDipoleContribution (water, dipolePotential);

		//printf ("%8.3e %8.3e\n", dipolePotential[0] * sfg_units::AU2WAVENUMBER, dipolePotential[1]*AU2WAVENUMBER);

		/***********************************************
		 * Adding in the dipole-dipole frequency shifts
		 ***********************************************/
		_OH1FreqShift += dipolePotential[0];
		_OH2FreqShift += dipolePotential[1];

		//printf ("% 10.3f\t% 10.3f\n", dipolePotential[0]*AU2WAVENUMBER, dipolePotential[1]*AU2WAVENUMBER);
		//}
#endif

		/*********************************************
		 * Final calculation of the OH bond frequency
		 *********************************************/
		// The two bond frequency values are now shifted from the gas phase vibrational value. (in cm-1)
		_w1 = UNCOUPLED_OH_FREQ + _OH1FreqShift;
		_w2 = UNCOUPLED_OH_FREQ + _OH2FreqShift;

		return;
}

// Here we calculate several values that play right off the equations of the paper - primarily Eq 7 and 8
/* such as the eigenfrequencies (normal modes) and the eigenvectors of the normal modes of the water. */
void SFGCalculator::WaterEigenSystem (morita::MoritaH2O& water) {

	// calculate the two frequency shifts of the OH bonds
	this->FreqShift (water);

	/********************
	 * Solving for symmetric and antisymmetric frequency of the normal modes
	 * ******************/

	double V12 = this->CouplingConstant (water);
	//double V12 = COUPLING_CONST;

	/* As per Dave's code in inter.f: */
	double wt = sqrt(_w1*_w1 + _w2*_w2 - 2.0*_w1*_w2 + 4.0*V12*V12);
	_ws = 0.5*(_w1+_w2 - wt) * 2.0*M_PI;
	_wa = 0.5*(_w1+_w2 + wt) * 2.0*M_PI;

	// let's set the symmetric frequency to be lower than the anti-symmetric
	/*
		 double temp = _ws;
		 if (_wa < _ws) {
		 _ws = _wa;
		 _wa = temp;
		 }
		 */

	/* Now again, as per Dave's solution: */
	_C1s = V12/(sqrt(V12*V12 + (_ws-_w1)*(_ws-_w1)));
	_C2s = (_ws-_w1)/(sqrt(V12*V12 + (_ws-_w1)*(_ws-_w1)));

	_C1a = V12/(sqrt(V12*V12 + (_wa-_w1)*(_wa-_w1)));
	_C2a = (_wa-_w1)/(sqrt(V12*V12 + (_wa-_w1)*(_wa-_w1)));

	/*
		 printf ("V12 = %12.4f\nwa = %12.4f ws = %12.4f\nw1 = %12.4f w2 = %12.4f\n", 
		 V12*AU2WAVENUMBER, _wa*AU2WAVENUMBER, _ws*AU2WAVENUMBER, _w1*AU2WAVENUMBER, _w2*AU2WAVENUMBER);
		 printf ("\nc1s = % 12.4f c1a = % 12.4f\nc2s = % 12.4f c2a = % 12.4f\n", _C1s, _C1a, _C2s, _C2a);
		 printf ("% 12.4f\n", (_w1-_wa)*_C1a + V12*_C2a);
		 exit(1);
		 */

	// before leaving, let's normalize the eigen vectors {C1x,C2x} to a magnitude of 1.0
	/*
		 double norm = sqrt(_C1s*_C1s + 1.0);
		 _C1s /= norm; _C2s /= norm;

		 norm = sqrt(_C1a*_C1a + 1.0);
		 _C1a /= norm; _C2a /= norm;

*/
	return;
}

/* Table 1 of the paper lists various ab initio properties of the water molecule for different axes/polarizations. Given an incoming polarization string (i.e. "SSP", "SPS", etc.) this function calculates the component of the dipole derivative for both OH bonds based on the rotation of the bond in the lab fixed-frame.

	 The polarization will be specified as S and P... we are dealing in molecular-frame coords, so we have to translate this to x, y, z. S means x and y because they are equivalent... however, you'll note that xx and yy polarization are not equivalent in the polarizability derivative - how do we address this discrepancy? The idea here is to calculate for both options. i.e. if the polarization desired is SPS, then we'll calculate both for xz and yz versions of the polarizability derivative, and then average over the two values... Let's see how well that works!

	 We're going to go through both OH bonds and treat each separately.
	 */
void SFGCalculator::PolarizabilityAndDipoleDerivs (morita::MoritaH2O& water, const int s1, const int s2, const int p) {

	// rotation matrix for moving from the water molecular frame to the lab frame.
	// let's find the rotation matrix for the water with which we're working
	// this matrix takes us from the frame of the first OH bond into the lab frame
	MatR _DCM (water.DCMToLabMorita(z));

	this->WaterEigenSystem (water);

	// Equation (14) from the paper adjusts the value of the dipole derivative term by scaling based on the frequency shift on the two OH bonds. Here we'll also calculate that scaling factor
	double scale1, scale2;

	// DSW applies the scaling factor to all components of the dipole deriv. vector
	// It seems like applying the scaling factor before or after rotation doesn't change much. It's just scaling, not changing direction, and so the constant multiplier gets carried throughout.
	// Also: Aren't these scaling factors supposed to be unitless? There shouldn't be any reason to convert back and forth from cm-1 or a.u. after calculating the factor as given in the M/H 2000 paper, equation (14)
	// here we need the frequency shift given as wavenumbers (cm-1).
	scale1 = (1.0 - 9.5138e-3 * _OH1FreqShift * sfg_units::AU2WAVENUMBER);
	scale2 = (1.0 - 9.5138e-3 * _OH2FreqShift * sfg_units::AU2WAVENUMBER);
	//printf ("% 10.3f% 10.3f\n", scale1, scale2);

	// New way of "pre-"rotating the alpha and mu tensors to avoid rotating the total Beta. a la DSW
	// here we'll set up the dipole and alpha derivative tensors, but rotated into the lab frame from where they were in the frame of the 1st OH bond

	// matrix rotation requires a unitary transformation to go into the lab frame
	//MatR rotAlpha1 = _DCM.transpose() * AlphaDeriv1 * _DCM;
	//MatR rotAlpha2 = _DCM.transpose() * AlphaDeriv2 * _DCM;
	MatR rotAlpha1 (_DCM * AlphaDeriv1 * _DCM.transpose());
	MatR rotAlpha2 (_DCM * AlphaDeriv2 * _DCM.transpose());

	// pull together all the combinations of the derivative terms for an average value
	_AlphaDerivS += _C1s * rotAlpha1(s1,s1) + _C2s * rotAlpha2(s1,s1);
	_AlphaDerivS += _C1s * rotAlpha1(s2,s2) + _C2s * rotAlpha2(s2,s2);
	// adding in the extra Alphaderiv cross-terms (x,z) & (z,x)...
	_AlphaDerivS += _C1s * rotAlpha1(s1,s2) + _C2s * rotAlpha2(s1,s2);
	_AlphaDerivS += _C1s * rotAlpha1(s2,s1) + _C2s * rotAlpha2(s2,s1);
	_AlphaDerivS /= 4.0;	 // and doing the requisite scaling

	_AlphaDerivA += _C1a * rotAlpha1(s1,s1) + _C2a * rotAlpha2(s1,s1);
	_AlphaDerivA += _C1a * rotAlpha1(s2,s2) + _C2a * rotAlpha2(s2,s2);
	_AlphaDerivA += _C1a * rotAlpha1(s1,s2) + _C2a * rotAlpha2(s1,s2);
	_AlphaDerivA += _C1a * rotAlpha1(s2,s1) + _C2a * rotAlpha2(s2,s1);
	_AlphaDerivA /= 4.0;

	// vector rotation into the lab frame is easy - just multiply by the rotation matrix
	VecR rotMu1 (_DCM * (MuDeriv1 * scale1));
	VecR rotMu2 (_DCM * (MuDeriv2 * scale2));

	_MuDerivS = _C1s * rotMu1[p] + _C2s * rotMu2[p];
	_MuDerivA = _C1a * rotMu1[p] + _C2a * rotMu2[p];

	//printf ("MuDeriv = \t%f\t%f\nAlphaDeriv = \t%f\t%f\n", _MuDerivS, _MuDerivA, _AlphaDerivS, _AlphaDerivA);
	return;
}


/* here we calculate the hyperpolarizability spectrum for a water molecule. In the course of this, two spectra will be calculated and averaged based on the two values of the eigenfrequencies of a water molecule */
std::vector< std::complex<double> >& SFGCalculator::Beta (morita::MoritaH2O& water, const int s1, const int s2, const int p) {
	//vector< complex<double> >& SFGCalculator::Beta (Water& water) {

	// for this step - the summation and rotation to the lab frame by use of the direction cosine matrix - see:
	// Shen - Phys. Rev. B, 59, 19 (1999), p. 12634
	// eq(7) of the paper lays out pretty simply that we use a direction cosine matrix to find all the interesting things that make up the hyperpolarizability.
	// l,m,n are the lab frame axes, and p,q,r are the molecular frame ones.

	/* Now that we have all the data established for this one water molecule (i.e. coefficients for symmetric and antisymmetric modes, dipole derivatives, and polarizability derivatives for all valid polarizations) we can put it all together and calculate spectra. Two spectra will come out for the sym and anti-sym cases, and we have to include all the polarization combinations. The beta that comes out of here is still in the molecular frame and needs to be rotated.
	*/

	// we've collected lots of information, so now let's construct some spectra!
	_Beta.clear();	// clear out the previous spectrum

	this->PolarizabilityAndDipoleDerivs (water, s1, s2, p);

	// let's calculate a pre-multiplier from equation (5)
	double premultSym = 0.5 / M / _ws * _AlphaDerivS * _MuDerivS;
	double premultAntiSym = 0.5 / M / _wa * _AlphaDerivA * _MuDerivA;

	for ( double wIR = START_FREQ; wIR <= END_FREQ; wIR += FREQ_STEP) {

		// a couple more handy values to compute
		double dFreq = _ws - wIR;
		double denominator = dFreq*dFreq + GAMMA_SQ;

		// and here we construct both real and imaginary parts of the beta spectrum
		double realPart, imagPart;

		// first calculate the symmetric part
		realPart = premultSym * dFreq / denominator;
		imagPart = premultSym * GAMMA / denominator;

		// and then recalculate to add in the anti-symmetric part as the response is the sum of both responses
		dFreq = _wa - wIR;
		denominator = dFreq*dFreq + GAMMA_SQ;

		realPart += premultAntiSym * dFreq/denominator;
		imagPart += premultAntiSym * GAMMA/denominator;

		_Beta.push_back (std::complex<double>(realPart, imagPart));
	}

	return (_Beta);
}

// calculates the dipole-dipole interaction potential for two dipoles separated by a distance R.
double SFGCalculator::DipolePotential (const VecR& muA, const VecR& muB, const VecR& R) {

	return 1.0/pow(R.Magnitude(), 3.0) * ((muA.dot(muB)) - 3.0/pow(R.Magnitude(), 2.0) * (muA.dot(R))*(R.dot(muB)));

}

// A scaling factor of 1/sqrt(N) is applied to the coupling constant to account for how OH intramolecular coupling is red-shifted when in solution. Thus, DSW introduced a factor of 1/sqrt(N), where N is the number of hydrogen bonds on the water molecule. 0 or 1 bond would leave the coupling constant unchanged.
double SFGCalculator::CouplingConstant (morita::MoritaH2O& water) const {

	double V12 = COUPLING_CONST;

	//	int coord = static_cast<int>(_matrix->WaterCoordination(&water));
	//int N = t.Graph().NumHBonds(&water);

	//if (N > 6) {
	//printf ("coordination == %d!\n", N);
	//exit(1);
	//}

	//if (N > 1) {
	//V12 /= sqrt((double)N);
	//}
	//printf ("coordination = %d% 10.3f\n", N, V12*AU2WAVENUMBER);

	return (V12);
}

#ifdef DIPOLE_DIPOLE
// This calculates the dipole-dipole contribution of a neighboring (source) water molecule to the (target) dipole potential
// Thus, the arguments are the target H through which H-bonding occurs, the pre-calculated target dipole derivative, and the dipole potential container that is being used to track the added contributions
void SFGCalculator::HDipoleDipoleContribution (
		const Atom * tH, 
		const VecR& t_mu, 
		const VecR& t_com, 
		const int oh_num,
		std::vector<double>& dipolePotential) 
{
	// Go through one of the target hydrogens and find their H-bonding partners
	Atom_ptr_vec hbonds = t.Graph().BondedAtoms(tH, hbond);
	for (Atom_it sO = hbonds.begin(); sO != hbonds.end(); sO++) {
		// find the source hydrogen, oxygen, center of mass, oh vector, etc.
		morita::MoritaH2O * sH2O = static_cast<morita::MoritaH2O *>((*sO)->ParentMolecule());
		VecR sOH1 (MDSystem::Distance(*sO, sH2O->GetAtom("H1")));		// in angstroms
		VecR sOH2 (MDSystem::Distance(*sO, sH2O->GetAtom("H2")));
		VecR sCOM1 ((*sO)->Position() + (sOH1.normalized() * (OH_COM_LENGTH/sfg_units::ANG2BOHR)));		// in angstroms
		VecR sCOM2 ((*sO)->Position() + (sOH2.normalized() * (OH_COM_LENGTH/sfg_units::ANG2BOHR)));

		// now we have to find the center of mass separation vectors, and also the dipole moment derivatives for each of the source OHs
		VecR sMu1 (sOH1.normalized() * MU_DERIV_LENGTH);		// atomic units
		VecR sMu2 (sOH2.normalized() * MU_DERIV_LENGTH);

		// The R vectors point from the source to the target
		VecR R1 (MDSystem::Distance(sCOM1, t_com));		// angstroms
		VecR R2 (MDSystem::Distance(sCOM2, t_com));

		// but we are dealing with atomic units, so let's rescale the R vectors to reflect this
		R1 *= sfg_units::ANG2BOHR;		// converting to atomic units
		R2 *= sfg_units::ANG2BOHR;

		dipolePotential[oh_num] += this->DipolePotential (t_mu, sMu1, R1);
		dipolePotential[oh_num] += this->DipolePotential (t_mu, sMu2, R2);
	}
	return;
}

void SFGCalculator::DipoleDipoleContribution (morita::MoritaH2O& water, std::vector<double>& dipolePotential) {
	/*********************************************************************************************
	 * Calculation of the frequency shift due to dipole-dipole interaction with neighboring waters
	 *********************************************************************************************/

	//int coord = static_cast<int>(graph->WaterCoordination(&water));
	//int N = graph->NumHBonds(&water);
	//if (N > 3) {
	// two OH bonds on the target water, so we will have two sets of shifts from dipole-dipole interactions

	// the dipole moment derivatives for each OH bond are calculated along the O->H vector, and have a magnitude of the same vector as in the MH paper.
	// For the dipole moment derivative, we can just take the one from the molecular frame and rotate it out to the lab frame
	// Originally the dipole moment derivative was assumed to be along the OH bond... as per DSW below:
	VecR oh1_unit (water.OH1()->normalized());
	VecR oh2_unit (water.OH2()->normalized());

	VecR mu1 (oh1_unit * MU_DERIV_LENGTH);	// these are in atomic units
	VecR mu2 (oh2_unit * MU_DERIV_LENGTH);
	// But perhaps another idea, if taking the actual dipole deriv. doesn't work, is to take the dot product of the actual one with the OH bond unit vector. (?)

	// find the center of mass locations for both target OH bonds
	VecR com1 (water.O()->Position() + (oh1_unit * OH_COM_LENGTH/sfg_units::ANG2BOHR));	// this is in Angstroms...
	VecR com2 (water.O()->Position() + (oh2_unit * OH_COM_LENGTH/sfg_units::ANG2BOHR));

	// Now we go through the calculation of each neighboring H-bonded water (the "source" waters, acting as the "source" of the dipole-dipole interactions) and find the contribution from each OH dipole. First we'll start with source OH's that are bound through the target oxygen. (I'll use the nomenclature of sH to mean source-hydrogen, and tO to mean target oxygen, etc.)

	Atom_ptr_vec hbonds = t.Graph().BondedAtoms(water.O(), hbond);	// all the H's H-bonding to the target O
	for (Atom_it sH = hbonds.begin(); sH != hbonds.end(); sH++) {

		// find the source hydrogen, oxygen, center of mass, oh vector, etc.
		morita::MoritaH2O * sH2O = static_cast<morita::MoritaH2O *>((*sH)->ParentMolecule());
		Atom * sO = sH2O->GetAtom("O");

		VecR sOH (MDSystem::Distance (sO, *sH));
		//VecR sOH (sO->Position().MinVector (sH->Position(), MDSystem::Dimensions()));	// source OH vector
		VecR sCOM (sO->Position() + (sOH.normalized() * (OH_COM_LENGTH/sfg_units::ANG2BOHR)));		// source OH center of mass in angstroms

		// now we have to find the center of mass separation vectors, and also the dipole moment derivatives for each of the source OHs
		VecR sMu (sOH.normalized() * MU_DERIV_LENGTH);		// in atomic units
		// The R vectors point from the source to the target
		VecR R1 (MDSystem::Distance(sCOM,com1));
		VecR R2 (MDSystem::Distance(sCOM,com2));

		// but we are dealing with atomic units, so let's rescale the R vectors to reflect this
		R1 *= sfg_units::ANG2BOHR;
		R2 *= sfg_units::ANG2BOHR;

		dipolePotential[0] += this->DipolePotential (mu1, sMu, R1);		// the result is in atomic units
		dipolePotential[1] += this->DipolePotential (mu2, sMu, R2);
	}

	HDipoleDipoleContribution (water.H1(), mu1, com1, 0, dipolePotential);
	HDipoleDipoleContribution (water.H2(), mu2, com2, 1, dipolePotential);
	//}

	return;
}

#endif
} // namespace morita
