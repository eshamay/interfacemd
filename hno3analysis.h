#ifndef HNO3ANALYSIS_H_
#define HNO3ANALYSIS_H_

#include "molecule.h"
#include "hno3.h"
#include "h2o.h"
#include "xyzsystem.h"

using namespace std;

const double R_OO = 3.5;			// These values are the maximums for defining H-bonds between the acceptor water and hno3.
const double R_OH = 2.45;			// This is the max value of an H-bond length (from rdf data)
const double A_HOO = 35.0;			// The angle of the H-O...O above which we don't have an H-bond.

const double R_OH_HNO3_mean = 1.016;  // avg OH distance on the NO3 according to Hynes... should be recalculate for each simulation.
const double R_OH_H3O_mean = 1.029;	  // same but for hydronium

// establish the value of the q_pt function as per Wang,Bianco,Hynes (2008). This is a function that "quantifies" the amount of dissociation of a contact ion pair. Given a nitric acid and an acceptor water molecule, this function returns a value that remains mostly between -1 and 1. A value of 0 means that the donor proton from the nitric acid is 1/2-way to dissociating to the acceptor water. A value of -1 is a complete nitric acid, and a value of 1 represents a pair of nitrate and hydronium. The contact-ion pair is what we end up with when the proton is somewhere in between the two such that the proton is "shared".
double q_pt (XYZSystem& sys, Molecule * hno3);

#endif
