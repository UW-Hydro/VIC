#include <stdio.h>
#include <stdlib.h>

void cshift(double *a, int width, int offset) 
{
   int i;
   double b;    
   
   b=*(a + width*offset);
   for (i = 0; i != width - 1; i++) {
      *(a + width*offset + i) = *(a + width*offset + i+1);
   }
   *(a + width*offset + i) = b;
   
    
}
