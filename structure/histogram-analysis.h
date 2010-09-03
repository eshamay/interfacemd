#ifndef HISTOGRAM_ANALYSIS_H_
#define HISTOGRAM_ANALYSIS_H_

#include "analysis.h"
#include "utility.h"

namespace md_analysis {

	template <typename T>
		// an analyzer for generating a single histogram (of doubles) of a system property
		class histogram_analyzer : public AnalysisSet< Analyzer<T> > {
			public:
				typedef Analyzer<T> system_t;

				histogram_analyzer (std::string desc, std::string fn) : AnalysisSet<system_t> (desc, fn) { }
				virtual ~histogram_analyzer() { }

				void DataOutput (system_t& t);
			protected:
				typedef std::vector<double> double_vec;
				double_vec values;
		};	// single histogram analyzer


	//template <>
		//histogram_analyzer<AmberSystem>::histogram_analyzer (std::string desc, std::string fn) : AnalysisSet<system_t> (desc,fn) { }
	//template <>
		//histogram_analyzer<XYZSystem>::histogram_analyzer (std::string desc, std::string fn) : AnalysisSet<system_t> (desc,fn) { }


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
