#ifndef __ANGLE_HISTOGRAM_H
#define __ANGLE_HISTOGRAM_H

#include "../utility.h"
#include "../analysis.h"

class AngleBinner : public unary_function<double,void>
{

  public:
  typedef Histogram1D<double> histo_t;
  void operator() (Water * wat)
  {
	/*
	wat->SetOrderAxes();
	//VecR v (wat->Z());				// set the bisector
	VecR v (wat->Y());
	histo(v < VecR(0.0,-1.0,0.0));
	*/
	wat->Print();
	wat->SetOrderAxes();
	VecR bis (wat->Z());				// set the bisector
	VecR norm (wat->Y());
	printf ("%d\t% 6.3f\t% 6.3f\n", wat->MolID(), bis < VecR(0.0,-1.0,0.0), norm < VecR(0.0,-1.0,0.0));
  }

  void Output (FILE * output)
  {
	rewind(output);
	for (double t = -1.0; t < 1.0; t += 0.01)
	{
	  fprintf(output, "% 8.3f %12d\n", t, histo.Population(t));
	}
  }

  private:
  static histo_t histo;

};


class AngleHistogram : public Analyzer {

  public:

	AngleHistogram (WaterSystemParams& wsp) :
	  Analyzer(wsp) { return; }

	AngleBinner binner;

	void LoadWaterSlice () {
	  LoadAll();
	  /* load only the waters */
	  FindWaters();
	  /* and only the ones within a certain slice of the interface */
	  SLICE_BY_POSITION(int_wats, Water *, 37.0, 38.5);
	}

	void Setup () { 
	  LoadWaterSlice();
	}

	void Analysis () {
	  LoadWaterSlice();
	  // bin each water's position and angle
	  FOR_EACH(int_wats, binner);
	  return;
	}

	// nothing at the moment
	void PostAnalysis () { return; }

	void DataOutput (const unsigned int timestep) {
	  if (!(timestep % (output_freq * 10)))
		binner.Output(output);
	  return;
	}
};

#endif
