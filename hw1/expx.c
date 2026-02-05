#include <petsc.h>
#include <math.h>

int main(int argc, char **argv) {
    PetscMPIInt    rank;
    PetscMPIInt    size; // number of processes
    PetscInt       i;
    PetscInt       n;
    PetscReal      x = 1.0, localsum, localproduct, globalsum;
    PetscInt       N = 1;
    PetscInt       termsPerProc, leftoverTerms;

    PetscCall(PetscInitialize(&argc,&argv,NULL,
        "Compute exp(x) in parallel with PETSc balanced for N-terms.\n\n"));
    PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
    PetscCall(MPI_Comm_size(PETSC_COMM_WORLD,&size));

    // read option
    PetscOptionsBegin(PETSC_COMM_WORLD,"","options for expx","");
    PetscCall(PetscOptionsReal("-x","input to exp(x) function",NULL,x,&x,NULL));
    PetscCall(PetscOptionsInt("-N","number of terms to use",NULL,N,&N,NULL));
    PetscOptionsEnd();

    termsPerProc = N / size;
    leftoverTerms = 0;
    localsum = 0.0;
    localproduct = 1.0;
    i = 1;

    // If there are too many processes, compute one term each for the first N ranks
    if (termsPerProc == 0) {
        termsPerProc = rank < N ? 1 : 0;
    } else if (rank == size - 1) {
        // Otherwise, if we're the last rank, take any leftover terms
        leftoverTerms = N % size;
    }

    // compute  x^n/n!
    for (n = 1 + rank * termsPerProc; n <= (rank + 1) * termsPerProc + leftoverTerms; n++) {
        PetscCall(PetscPrintf(PETSC_COMM_SELF,
            "rank %d computing term n=%d\n",rank,n));
        for (; i < n; i++) { // Don't reset to 1, save work by starting from previous product
            localproduct *= x/i;
        }

        localsum += localproduct;
    }

    // sum the contributions over all processes
    PetscCall(MPI_Allreduce(&localsum,&globalsum,1,MPIU_REAL,MPIU_SUM,
        PETSC_COMM_WORLD));

    // output estimate and report on work from each process
    if (rank == 0) { // Ensure process 0 does the printing
        unsigned long error = fabs(exp(x) - globalsum)/PETSC_MACHINE_EPSILON;

        PetscCall(PetscPrintf(PETSC_COMM_SELF,
            "exp(%17.15f) is about %17.15f (error = %ld*EPS)\n",x,globalsum,error));
    }

    PetscCall(PetscFinalize());
    return 0;
}
