void swap(double *a, double *b){
double temp = *a;
*a = *b;
*b = temp;
}

int partition(double *A, int left, int right){
   double pivot = A[right];
   int i = left;
   int x;
   for (x = left; x < right; x++) {
      if (A[x] < pivot){
         swap(&A[i], &A[x]);
         i++;
      }
   }
   swap(&A[i], &A[right]);
   return i;
}


double* quickselect(double *A, int left, int right, int k){
//p is position of pivot in the partitioned array
   int p = partition(A, left, right);
//k equals pivot got lucky
   if (p == k-1){
      double *temp = malloc((k)*sizeof(double));
      for(int i=left; i<=k-1; ++i) {
         temp[i]=A[i];
      }
      return temp;
   }
//k less than pivot
   else if (k - 1 < p){
      return quickselect(A, left, p - 1, k);
   }
//k greater than pivot
   else{
      return quickselect(A, p + 1, right, k);
   }
}
