BEGIN 		{i = 0}
$NR == 1 	{print}
$6 != i 	{print "TER";
			 i = $6}
$6 == i 	{print}
