BEGIN { 
  FS = ",";
  i = 0; 
  while (getline var < "cosby.1984.csv") {
    n=split(var,var2,","); 
    for (j = 1; j <= n; j++) { 
      array[i,j] = var2[j];
    } 
    i++;
  }
  nclasses = i-1;
}
{
  if ($3 == "unknown texture")
    print "-1 -1 -1 -1 -1";
  else {
    id = 0;
    Ks = 0;
    b = 0;
    Ts = 0;
    Fc = 0;
    Wp = 0;
    n = split($3,var,"/");
    for (i = 1; i <= n; i++) {
      for (j = 1; j <= nclasses; j++) {
	if (var[i] == array[j,1]) {
	  id += j;
	  Ks += array[j,3];
	  b += array[j,4];
	  Ts += array[j,5];
	  Fc += array[j,7];
	  Wp += array[j,8];
	  break;
	}
      }
    }
    id /= n;
    Ks /= n;
    b /= n;
    Ts /= n;
    Fc /= n;
    Wp /= n;
    print Ks, Ts, Fc, Wp, b;
  }
}
