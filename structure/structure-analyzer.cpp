#include "structure-analyzer.h"

int main () {

	printf ("Choose the system type to analyze from the list below\n\n");
	typedef enum {AMBER, XYZ} system_type;
	printf ("%d) Amber System\n%d) XYZ System\n\n", AMBER, XYZ);
	printf ("analysis choice --> ");
	int choice;
	std::cin >> choice;

	if (choice == AMBER)
		md_analysis::StructureAnalyzer<AmberSystem> sa;
	else if (choice == XYZ)
		md_analysis::StructureAnalyzer<XYZSystem> sa;

	return 0;
}


