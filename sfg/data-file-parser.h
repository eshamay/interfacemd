#ifndef DATA_FILE_PARSER_H_
#define DATA_FILE_PARSER_H_

// check out http://mybyteofcode.blogspot.com/2010/02/parse-csv-file-with-boost-tokenizer-in.html
// info on the boost tokenizer
#include "matrixr.h"
#include "utility.h"
#include <map>
#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <functional>
#include <iterator>     // ostream_operator

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map.hpp>


namespace datafile_parsers {

	using namespace boost;

	const double multiplier = 1000.0;

	class DataFile {
		public:
			DataFile (const std::string filename) : _filename(filename), _input(_filename.c_str()), _numLines(0) {
				if (!_input.is_open()) {
					std::cout << "Couldn't open the file " << _filename << std::endl;
					exit(1);
				}
			}

			virtual ~DataFile () {
				if (_input.is_open())
					_input.close();
			}

		protected:
			std::string _filename;
			std::ifstream _input;
			long int _numLines;
	};




	//class DataParser : public std::unary_function <std::vector<
	class DelimitedDataFile : public DataFile {
		public:
			DelimitedDataFile (const std::string filename, const char* delimiter) : DataFile(filename), _delimiter(delimiter) { 
			}

			virtual void ParseFile ();	// parses through an entire file
			virtual void ParseRow (const std::vector<double>&) = 0;		// takes care of processing each row in the file
			virtual void ParseData () = 0;		// manipulates the total data collected from each row

		protected:
			typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;

			typedef boost::char_separator<char> Delimiter;
			Delimiter _delimiter;
	};	// Data file parser




	template <class T>
		class MultiKeyValueGraph {
			public:
				MultiKeyValueGraph () { 
					_root = add_vertex(_graph);
				}

				void Insert (const double r1, const double r2, const double theta, const T& t);
				T& Value (const double r1, const double r2, const double theta);

			protected:
				struct VertexProperties {
					int key;
					T value;
				};

				typedef adjacency_list<listS, listS, directedS, VertexProperties> Graph;
				typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
				typedef typename graph_traits<Graph>::vertex_iterator Vertex_it;
				typedef std::pair<Vertex_it, Vertex_it> Vertex_pair;
				typedef typename graph_traits<Graph>::edge_descriptor Edge;
				typedef typename graph_traits<Graph>::edge_iterator Edge_it;
				typedef typename graph_traits<Graph>::adjacency_iterator Adj_it;

				// generic property maps
				template <class U, class Property_U> struct PropertyMap { typedef typename property_map<Graph, U Property_U::*>::type Type; };

				static typename PropertyMap<T,VertexProperties>::Type 		v_value;
				static typename PropertyMap<int,VertexProperties>::Type			v_key;

				Vertex AdjacentVertexAddValue (Vertex& start, const double key);
				Vertex FindVertex (Vertex& start, const double key);

				static Graph _graph;
				Vertex _root;
		};


	template <class T>
		MultiKeyValueGraph<T>::Graph MultiKeyValueGraph<T>::_graph(0);

	template <class T>
		typename MultiKeyValueGraph<T>::PropertyMap<T,MultiKeyValueGraph<T>::VertexProperties>::Type MultiKeyValueGraph<T>::v_value = get(&VertexProperties::value, _graph);

	template <class T>
		typename MultiKeyValueGraph<T>::PropertyMap<int,MultiKeyValueGraph<T>::VertexProperties>::Type MultiKeyValueGraph<T>::v_key = get(&VertexProperties::key, _graph);

	template <class T>
		MultiKeyValueGraph<T>::Vertex MultiKeyValueGraph<T>::AdjacentVertexAddValue (Vertex& start, const double key) {

			//printf ("AdjacentVertexAddValue: adding %5.2f to %5.2f\n", key, v_key[start]);
			Vertex vd;
			Edge e;
			bool added;
			Adj_it vi, vi_end, next;
			tie(vi, vi_end) = adjacent_vertices(start, _graph);
			int val = int(key*multiplier);
			for (next = vi; vi != vi_end; vi = next) {
				next++;

				if (v_key[*vi] == val) {
					//printf ("found %5.2f already connected\n", v_key[*vi]);
					return *vi;
				}
			}
			if (vi == vi_end) {
				// if the vertex doesn't already exist, then add it and tie it to the parent
				vd = add_vertex(_graph);
				v_key[vd] = int(key*multiplier);
				//printf ("now tying the two keys together: %5.2f --> %5.2f\n", v_key[start], v_key[vd]);
				tie(e,added) = add_edge(start,vd,_graph);
				return vd;
			}
			else {
				//printf ("found key: %5.2f\n", v_key[*vi]);
				return *vi;
			}
		}

	template <class T>
		MultiKeyValueGraph<T>::Vertex MultiKeyValueGraph<T>::FindVertex (Vertex& start, const double key) {
			//printf ("FindVertex: searching for %5.2f starting at %5.2f\n", key, v_key[start]);
			Adj_it vi, vi_end, next;
			tie(vi, vi_end) = adjacent_vertices(start, _graph);
			int val = int(key*multiplier);
			//printf ("val=%d\n", val);
			for (next = vi; vi != vi_end; vi = next) {
				next++;
				//printf ("found %5.2f\n", v_key[*vi]);

				if (v_key[*vi] == val) {
					return *vi;
				}
			}
			if (vi == vi_end) { 
				std::cout << "couldn't find the vertex for key: " << key << std::endl;
				exit(1);
			}
			return *vi;
		}

	template <class T>
		T& MultiKeyValueGraph<T>::Value (const double r1, const double r2, const double theta) {
			Vertex vd = FindVertex (_root, r1);
			vd = FindVertex (vd, r2);
			vd = FindVertex (vd, theta);
			return v_value[vd];
		}


	template <class T>
		void MultiKeyValueGraph<T>::Insert (const double r1, const double r2, const double theta, const T& mat) {
			//printf ("Insert: adding %5.2f %5.2f %6.2f\n", r1, r2, theta);
			Vertex vd = AdjacentVertexAddValue (_root, r1);
			vd = AdjacentVertexAddValue (vd, r2);
			vd = AdjacentVertexAddValue (vd, theta);
			v_value[vd] = mat;
		}



	template <class T>
		class ElectrostaticMomentFile : public DelimitedDataFile {
			public:
				ElectrostaticMomentFile (const std::string filename)
					: DelimitedDataFile(filename, ";"), r_min(0.93), r_max(1.13), dr(0.02), theta_min(90.0), theta_max(120.0), dtheta(3.0)
					{ }

				T& Value (const double r1, const double r2, const double theta);

			protected:
				MultiKeyValueGraph<T> _graph;
				double r_min;	// bondlength minimum key key
				double r_max;
				double dr;	// resolution of data points along the bondlengths
				double theta_min;
				double theta_max;
				double dtheta;

				virtual void ParseRow (const std::vector<double>&) = 0;
				virtual void ParseData () { }
				//void Polarizability (double r1, const double r2, const double theta);
				virtual double MatchValue (const double key, const double resolution, const double min, const double max);
		};




	class DipoleDataFile : public ElectrostaticMomentFile<VecR> {
		public:
			DipoleDataFile (const std::string filename) : ElectrostaticMomentFile<VecR> (filename) { this->ParseFile(); }
		protected:
			void ParseRow (const std::vector<double>&);
	};




	class PolarizabilityDataFile : public ElectrostaticMomentFile<MatR> {
		public:
			PolarizabilityDataFile (const std::string filename) : ElectrostaticMomentFile<MatR> (filename) { this->ParseFile(); }
			//! The Value method returns the polarizability for a given water geometry (r1, r2, theta). The returned matrix is in the reference frame designated in Morita/Hynes 2000 - the Z-axis is along one of the bonds (r1) and the other OH bond is in the positive x-axis direction, and the molecule is in the x-z plane.
		protected:
			void ParseRow (const std::vector<double>&);
	};	// polarizability data file



	template <typename T>
		class StringToNumber : public std::unary_function<std::string,T> {
			public:
				T operator() (const std::string& str) const {
					return boost::lexical_cast<T>(str);
				}
		};



	void DelimitedDataFile::ParseFile () {

		std::vector< double > vec;
		std::string line;

		while (getline(_input,line)) {
			Tokenizer tok(line,_delimiter);
			std::transform(tok.begin(), tok.end(), std::back_inserter(vec), StringToNumber<double>());

			this->ParseRow(vec);
			_numLines++;

			vec.clear();
		}

		this->ParseData();
	}


	void DipoleDataFile::ParseRow (const std::vector<double>& vec) {

		std::vector<double>::const_iterator it = vec.begin();
		double r1 = *it++;
		double r2 = *it++;
		double theta = *it++;

		double m[3];
		std::copy(it, vec.end(), m);
		VecR v (m);	// the matrix extracted from each row

		//printf ("ParseRow: adding %5.2f %5.2f %6.2f\n", r1, r2, theta);
		_graph.Insert (r1, r2, theta, v);
	}	// parse row


	void PolarizabilityDataFile::ParseRow (const std::vector<double>& vec) {

		double m[9];
		std::vector<double>::const_iterator it = vec.begin();
		double r1 = *it++;
		double r2 = *it++;
		double theta = *it++;
		std::copy(it, vec.end(), m);
		MatR mat (m);	// the matrix extracted from each row

		//printf ("ParseRow: adding %5.2f %5.2f %6.2f\n", r1, r2, theta);
		_graph.Insert (r1, r2, theta, mat);

	}	// parse row


	//! Take an incoming key and fix it such that it matches one of the keys in the dataset
	template <class T>
	double ElectrostaticMomentFile<T>::MatchValue (const double key, const double resolution, const double min, const double max) {

		int step = int((key*multiplier-min*multiplier)/multiplier/resolution);
		int floor = resolution*multiplier*step+min*multiplier;
		floor = (floor > int(max*multiplier)) ? int(max*multiplier) : floor;	// a bounds check
		floor = (floor < int(min*multiplier)) ? int(min*multiplier) : floor;	// a bounds check

		return floor/multiplier;
	}


	template <class T>
	T& ElectrostaticMomentFile<T>::Value (const double r1, const double r2, const double theta) {
		// first cast the numbers to the proper keys
		double r1_floor = MatchValue(r1, dr, r_min, r_max);
		double r2_floor = MatchValue(r2, dr, r_min, r_max);
		double theta_floor = MatchValue(theta, dtheta, theta_min, theta_max);

		//printf ("%f %f %f\n", r1, r2, theta);
		//printf ("%f %f %f\n", r1_floor, r2_floor, theta_floor);
		fflush(stdout);
		return _graph.Value(r1_floor, r2_floor, theta_floor);
	}



}	 // namespace datafile parsers

#endif
