#ifndef SFGUNITS_H_
#define SFGUNITS_H_

namespace sfg_units {


  /* Note: I think morita-hynes use atomic units (bohrs, hartrees, and the like)  so here are various units and conversion factors */

  // in atomic units, this is the mass of a proton. Mass of electron = 1.0
  const  double MPROT	=	1836.153;				
  // mass of an oxygen atom in atomic units
  const  double MOXY	=	(MPROT*8.0+8.0);		
  // and mass of a hydrogen
  const  double MHYD	=	(MPROT*1.0+1.0);		

  // some very useful conversion factors
  // angstroms to bohr radii
  const double ANG2BOHR			=	1.8897161646320724;					

  // from hartree to kcal/mol
  const double HARTREE2KCALPMOL	=	627.509;						
  // convert amber forces (kcal/mol/A) into atomic force units
  const double AMBER2ATOMIC		=	1.0/HARTREE2KCALPMOL/ANG2BOHR;	

  // convert from Hz to wavenumbers (cm-1)  (this is 1/c)
  const double HZ2WAVENUMBER	=	3.335641e-11;				
  // convert Hz to atomic units of frequency
  const double HZ2AU			=	2.418884324306202e-17*2.0*M_PI;			
  // convert from frequencies in atomic units to cm-1 (note: **not angular frequencies!** For that we need to fix the factor of 2*Pi)
  const double AU2WAVENUMBER		=	HZ2WAVENUMBER/HZ2AU;			

} // namespace sfg_units

#endif
