#ifndef HISTOGRAM_ANALYSIS_H_
#define HISTOGRAM_ANALYSIS_H_

#include "analysis.h"
#include "utility.h"

namespace md_analysis {

	template <typename T, typename U=double>
		class histogram_analysis : public AnalysisSet< Analyzer<T> > {
			public:
				typedef Analyzer<T> system_t;

				histogram_analysis (
						std::string desc, std::string fn,
						const U min, const U max, const U resolution)
					: 
						AnalysisSet<system_t> (desc, fn),
						histogram (min, max, resolution) { }

				virtual ~histogram_analysis() { }

			protected:
				typedef histogram_utilities::Histogram1D<U>	histogram_t;
				histogram_t histogram;
		};



	template <typename T>
		class double_histogram_analysis : public histogram_analysis<T> {
			public:
				typedef typename histogram_analysis<T>::system_t system_t;

				double_histogram_analysis (
						std::string desc, std::string fn,
						const double min, const double max, 
						const double min_2, const double max_2, const int number_of_bins)
					: 
						histogram_analysis<T> (desc, fn, min, max, (max-min)/double(number_of_bins)),
						histogram_2 (min_2, max_2, (max_2-min_2)/double(number_of_bins)) { }

				virtual ~double_histogram_analysis() { }
				virtual void DataOutput (system_t&);

			protected:
				typedef typename histogram_analysis<T>::histogram_t histogram_t;
				histogram_t histogram_2;

		};



	template <typename T>
		void double_histogram_analysis<T>::DataOutput(system_t& t) {

			rewind(t.Output());

			for (double val_1 = this->histogram.Min(), val_2 = this->histogram_2.Min(); val_1 < this->histogram.Max(); 
					val_1 += this->histogram.Resolution(), val_2 += this->histogram_2.Resolution()) {
				fprintf (t.Output(), "% 8.3f % 8.3f % 8.3f % 8.3f\n", val_1, double(this->histogram.Population(val_1))/double(t.Timestep()), val_2, double(this->histogram_2.Population(val_2))/double(t.Timestep()));
			}

			fflush(t.Output());
		}



	template <typename T>
		// an analyzer for generating a single histogram (of doubles) of a system property
		class histogram_analyzer : public AnalysisSet< Analyzer<T> > {
			public:
				typedef Analyzer<T> system_t;

				histogram_analyzer (std::string desc, std::string fn) : AnalysisSet<system_t> (desc, fn) { }
				virtual ~histogram_analyzer() { }

				virtual void DataOutput (system_t& t);
			protected:
				typedef std::vector<double> double_vec;
				double_vec values;
		};	// single histogram analyzer


	template <typename T>
		// an analyzer for generating histograms (of doubles) of two system properties
		class double_histogram_analyzer : public histogram_analyzer<T> {
			public:
				typedef typename histogram_analyzer<T>::system_t system_t;

				double_histogram_analyzer (std::string desc, std::string fn) : histogram_analyzer<T> (desc,fn) { }
				virtual ~double_histogram_analyzer() { }
				void DataOutput (system_t& t);
			protected:
				typedef typename histogram_analyzer<T>::double_vec double_vec;
				double_vec second_values;
		};	// double histogram analyzer




	template <typename T>
		void histogram_analyzer<T>::DataOutput (system_t& t) {

			rewind(t.Output());

			histogram_utilities::histogram_t first_histo = histogram_utilities::Histogram (values.begin(), values.end(), 100);
			int first_histo_max = histogram_utilities::MaxPopulation (first_histo.begin(), first_histo.end());

			for (unsigned i = 0; i < first_histo.size(); i++) {
				fprintf (t.Output(), "% 8.4f % 8.4f\n",
						first_histo[i].first, (double)first_histo[i].second/first_histo_max);
			}
			return;
		}


	template <typename T>
		void double_histogram_analyzer<T>::DataOutput (system_t& t) {

			rewind(t.Output());

			histogram_utilities::histogram_t first_histo = histogram_utilities::Histogram (this->values.begin(), this->values.end(), 200);
			int first_histo_max = histogram_utilities::MaxPopulation (first_histo.begin(), first_histo.end());

			histogram_utilities::histogram_t second_histo = histogram_utilities::Histogram (this->second_values.begin(), this->second_values.end(), 200);
			int second_histo_max = histogram_utilities::MaxPopulation (second_histo.begin(), second_histo.end());

			for (unsigned i = 0; i < first_histo.size(); i++) {
				fprintf (t.Output(), "% 8.4f % 8.4f % 8.4f % 8.4f\n", 
						first_histo[i].first, (double)first_histo[i].second/first_histo_max, 
						second_histo[i].first, (double)second_histo[i].second/second_histo_max);
			}
			return;
		}

}	// namespace md_analysis

#endif
