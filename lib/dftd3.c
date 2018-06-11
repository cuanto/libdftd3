void dftd3(const int natoms, const double *coords, const int *itype, 
           const char *func, const int version, const int tz,
           double *edisp, double *grad){
  wrapper_(&natoms, coords, itype, func, &version, &tz, edisp, grad);
}
