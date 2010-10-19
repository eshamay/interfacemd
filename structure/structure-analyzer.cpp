#include "structure-analyzer.h"

int main (int argc, char**argv) {

	typedef enum {AMBER, XYZ} system_type;

	if (argc < 2) {
		printf ("Run this program using the following syntax:\n");
		printf ("structure-analyzer <system-type>\n\n");
		printf ("%d) Amber System\n%d) XYZ System\n\n", AMBER, XYZ);
		exit(1);
	}

	//printf ("analysis choice --> ");
	int system_choice = atoi(argv[1]);
	int analysis_choice = -1;
	if (argc == 3)
		analysis_choice = atoi(argv[2]);

	if (system_choice == AMBER)
		md_analysis::StructureAnalyzer<AmberSystem> sa(analysis_choice);
	else if (system_choice == XYZ)
		md_analysis::StructureAnalyzer<XYZSystem> sa(analysis_choice);

	return 0;
}


