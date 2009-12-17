{
	if (NF < 3) {
		print
	}
	else
	{
		printf ("%-6s%5d %-4s %3s  %4d    %8.3f%8.3f%8.3f\n", $1, $2, $3, $4, $5, $6, $7, $8);
	}
}
