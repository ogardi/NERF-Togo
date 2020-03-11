#include <R.h> 

void phcf_filter(int *nrows, int *ncols, int *x){
  
  int nr = nrows[0];
  
  int nc = ncols[0];
  
  int i;
  
  // first loop: keep all defor pixels in a square of 4 to 1, all others to 5
  for(i=nc+1; i<nr*nc-nc-1; i++) {
    // not considering first and last rows / columns
    if(i%nc!=0 && i%nc!=1) {
      if (x[i]==1) {
        if (x[i-nc]==1 && x[i-nc+1]==1 && x[i+1]==1 || x[i+1]==1 && x[i+nc+1]==1 && x[i+nc]==1 || x[i-nc-1]==1 && x[i-nc]==1 && x[i-1]==1 ||x[i-1]==1 && x[i+nc]==1 && x[i+nc-1]==1) {
          x[i] = 1;
        } else {
          x[i] = 5;
        }
      }
    }
  }
  
  // second loop, set all defor-pixels of type 5 adjacent to a type 1 pixel to 6
  for(i=nc+1; i<nr*nc-nc-1; i++) {
    if(i%nc!=0 && i%nc!=1) {
      if (x[i]==5) {
        if (x[i-nc]==1 || x[i+nc]==1 || x[i+1]==1 || x[i-1]==1) {
          x[i] = 6;
        }
      }
    }
  }
  
  // third loop, set all defor-pixels of type 6 back to 1 
  // set all remainingg defor-pixels of type 5 either to 7 (forest) or 0 (non-forest)
  for(i=nc+1; i<nr*nc-nc-1; i++) {
    if(i%nc!=0 && i%nc!=1) {
      if (x[i]==6) {
        x[i] = 1;
      }
      if (x[i]==5) {
        if (x[i-nc]==2 || x[i+nc]==2 || x[i+1]==2 || x[i-1]==2) {
          x[i] = 7;
        } else {
          x[i] = 0;
        }
      }
    }
  }
  
  // fourth loop, set all forest pixels of type 7 back to 2 
  for(i=nc+1; i<nr*nc-nc-1; i++) {
    if(i%nc!=0 && i%nc!=1) {
      if (x[i]==7) {
        x[i] = 2;
      }
    }
  }
  
}


