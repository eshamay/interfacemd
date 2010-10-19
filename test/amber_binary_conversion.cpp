#include "amber_binary_conversion.h"

int main (int argc, char **argv) {

	if (argc < 3) {
		printf ("no timestep specified!\n");
		exit(1);
	}
	AmberSystem sys ("prmtop", "mdcrd");

	FILE * pfile;
	pfile = fopen("test.bin", "wb");

	int timesteps = atoi(argv[1]);
	float positions[sys.size() * 3];
	VecR pos;

	for (int i = 0; i < timesteps; i++) {
		for (int a = 0; a < sys.size(); a++) {
			pos = sys[a]->Position();
			positions[3*a] = float(pos[x]);
			positions[3*a+1] = float(pos[y]);
			positions[3*a+2] = float(pos[z]);
		}

		fwrite (&positions, sizeof(float), sys.size()*3, pfile);
	}

	fclose(pfile);

	timesteps = atoi(argv[2]);
	pfile = fopen("test.bin", "rb");
	for (int i = 0; i < timesteps; i++) {
		fread (&positions, sizeof(float), sys.size()*3, pfile);

		for (int a = 0; a < sys.size(); a++) {
			printf ("% 12.7f % 12.7f % 12.7f\n", positions[3*a], positions[3*a+1], positions[3*a+2]);
		}
	}

	fclose(pfile);
	return 0;
}
