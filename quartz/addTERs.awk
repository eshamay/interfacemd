BEGIN { i = i }

{
	if (NF < 3) { print }
	else if ($5 == i) { print }
	else if ($5 != i)
	{
		print;
		print "TER";
		i = $5
	}
}

