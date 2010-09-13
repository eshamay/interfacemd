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


	class MultiKeyValueGraph {
		public:
			MultiKeyValueGraph () { 
				_root = add_vertex(_graph);
			}

			void Insert (const double r1, const double r2, const double theta, const MatR& mat);
			MatR& Matrix (const double r1, const double r2, const double theta);

		protected:
			struct VertexProperties {
				int value;
				MatR matrix;
			};

			typedef adjacency_list<listS, listS, directedS, VertexProperties> Graph;
			typedef graph_traits<Graph>::vertex_descriptor Vertex;
			typedef graph_traits<Graph>::vertex_iterator Vertex_it;
			typedef std::pair<Vertex_it, Vertex_it> Vertex_pair;
			typedef graph_traits<Graph>::edge_descriptor Edge;
			typedef graph_traits<Graph>::edge_iterator Edge_it;
			typedef graph_traits<Graph>::adjacency_iterator Adj_it;

			// generic property maps
			template <class T, class Property_T> struct PropertyMap { 
				typedef typename property_map<Graph, T Property_T::*>::type Type; };


			static PropertyMap<MatR,VertexProperties>::Type 		v_matrix;
			static PropertyMap<int,VertexProperties>::Type 		v_value;

			Vertex AdjacentVertexAddValue (Vertex& start, const double value);
			Vertex FindVertex (Vertex& start, const double value);

			static Graph _graph;
			Vertex _root;
	};


	MultiKeyValueGraph::Graph MultiKeyValueGraph::_graph(0);

	MultiKeyValueGraph::PropertyMap<MatR,MultiKeyValueGraph::VertexProperties>::Type MultiKeyValueGraph::v_matrix = get(&VertexProperties::matrix, _graph);

	MultiKeyValueGraph::PropertyMap<int,MultiKeyValueGraph::VertexProperties>::Type MultiKeyValueGraph::v_value = get(&VertexProperties::value, _graph);


	MultiKeyValueGraph::Vertex MultiKeyValueGraph::AdjacentVertexAddValue (Vertex& start, const double value) {

		//printf ("AdjacentVertexAddValue: adding %5.2f to %5.2f\n", value, v_value[start]);
		Vertex vd;
		Edge e;
		bool added;
		Adj_it vi, vi_end, next;
		tie(vi, vi_end) = adjacent_vertices(start, _graph);
		int val = int(value*100.0);
		for (next = vi; vi != vi_end; vi = next) {
			next++;

			if (v_value[*vi] == val) {
				//printf ("found %5.2f already connected\n", v_value[*vi]);
				return *vi;
			}
		}
		if (vi == vi_end) {
			// if the vertex doesn't already exist, then add it and tie it to the parent
			vd = add_vertex(_graph);
			v_value[vd] = int(value*100.0);
			//printf ("now tying the two values together: %5.2f --> %5.2f\n", v_value[start], v_value[vd]);
			tie(e,added) = add_edge(start,vd,_graph);
			return vd;
		}
		else {
			//printf ("found value: %5.2f\n", v_value[*vi]);
			return *vi;
		}
	}

	MultiKeyValueGraph::Vertex MultiKeyValueGraph::FindVertex (Vertex& start, const double value) {
		//printf ("FindVertex: searching for %5.2f starting at %5.2f\n", value, v_value[start]);
		Adj_it vi, vi_end, next;
		tie(vi, vi_end) = adjacent_vertices(start, _graph);
		int val = int(value*100.0);
		//printf ("val=%d\n", val);
		for (next = vi; vi != vi_end; vi = next) {
			next++;
			//printf ("found %5.2f\n", v_value[*vi]);

			if (v_value[*vi] == val) {
				return *vi;
			}
		}
		if (vi == vi_end) { 
			std::cout << "couldn't find the vertex for value: " << value << std::endl;
			exit(1);
		}
		return *vi;
	}

	MatR& MultiKeyValueGraph::Matrix (const double r1, const double r2, const double theta) {
		Vertex vd = FindVertex (_root, r1);
		vd = FindVertex (vd, r2);
		vd = FindVertex (vd, theta);
		return v_matrix[vd];
	}


	void MultiKeyValueGraph::Insert (const double r1, const double r2, const double theta, const MatR& mat) {
		//printf ("Insert: adding %5.2f %5.2f %6.2f\n", r1, r2, theta);
		Vertex vd = AdjacentVertexAddValue (_root, r1);
		vd = AdjacentVertexAddValue (vd, r2);
		vd = AdjacentVertexAddValue (vd, theta);
		v_matrix[vd] = mat;
	}



	class PolarizabilityDataFile : public DelimitedDataFile {
		public:
			PolarizabilityDataFile (const std::string filename) 
				: 
					DelimitedDataFile(filename, ";"), r_min(0.93), r_max(1.13), dr(0.02), theta_min(90.0), theta_max(120.0), dtheta(3.0) 
			{
				this->ParseFile();
			}


			//! The Matrix method returns the polarizability for a given water geometry (r1, r2, theta). The returned matrix is in the reference frame designated in Morita/Hynes 2000 - the Z-axis is along one of the bonds (r1) and the other OH bond is in the positive x-axis direction, and the molecule is in the x-z plane.
			MatR& Matrix (const double r1, const double r2, const double theta);

		protected:
			MultiKeyValueGraph _graph;
			double r_min;	// bondlength minimum key value
			double r_max;
			double dr;	// resolution of data points along the bondlengths
			double theta_min;
			double theta_max;
			double dtheta;

			void ParseRow (const std::vector<double>&);
			void ParseData () { }
			//void Polarizability (double r1, const double r2, const double theta);
			double MatchValue (const double value, const double resolution, const double min, const double max);
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


	//! Take an incoming value and fix it such that it matches one of the keys in the dataset
	double PolarizabilityDataFile::MatchValue (const double value, const double resolution, const double min, const double max) {

		int step = int((value*100.0-min*100.0)/100.0/resolution);
		int floor = resolution*100.0*step+min*100.0;
		floor = (floor > int(max*100.0)) ? int(max*100.0) : floor;	// a bounds check
		floor = (floor < int(min*100.0)) ? int(min*100.0) : floor;	// a bounds check
		
		return floor/100.0;
	}


	MatR& PolarizabilityDataFile::Matrix (const double r1, const double r2, const double theta) {
		// first cast the numbers to the proper values
		double r1_floor = MatchValue(r1, dr, r_min, r_max);
		double r2_floor = MatchValue(r2, dr, r_min, r_max);
		double theta_floor = MatchValue(theta, dtheta, theta_min, theta_max);

		//printf ("%f %f %f\n", r1, r2, theta);
		//printf ("%f %f %f\n", r1_floor, r2_floor, theta_floor);
		fflush(stdout);
		return _graph.Matrix(r1_floor, r2_floor, theta_floor);
	}


	/*
	void PolarizabilityDataFile::Polarizability (double r1, const double r2, const double theta) {
		printf ("Polarizability: searching for %8.4f %8.4f %8.4f\n", r1, r2, theta);
		MatR mat = this->Matrix(r1,r2,theta);
		mat.Print();
	}
	*/

}	 // namespace datafile parsers

#endif
