#include <petsc.h>
#include <math.h>

int main(int argc, char **argv) {
    PetscMPIInt    rank;
    PetscMPIInt    size; // number of processes
    PetscInt       i, k;
    PetscInt       n;
    PetscReal      x = 1.0, localsum, localprod, globalsum, multiplier, prevMultiplier;
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
    localsum = 1.0;
    localprod = 1.0;
    k = rank * termsPerProc;
    
    for (i = 1; i < termsPerProc; i++) { // i = local term index
        k = rank * termsPerProc + i; // k = global term index

        localprod *= x / k;
        localsum += localprod;
        
        // For debugging only
        // PetscCall(PetscPrintf(PETSC_COMM_SELF,
        //     "rank %d: term %d: x / %d\n",rank,k,k));
    }

    multiplier = localprod * x / (k + 1); // Compute one more term for stitching

    // For debugging only
    // PetscCall(PetscPrintf(PETSC_COMM_SELF,
    //     "rank %d: multiplier for stitching last term: x / %d\n",rank,k + 1));

    // Stitching process
    // Note: it's ok to wait for message here since the heavy work is done
    // so we're not blocking any expensive computation
    if (rank > 0) {
        // Receive multiplier from previous rank
        PetscCall(MPI_Recv(&prevMultiplier, 1, MPIU_REAL, rank - 1, 0,
            PETSC_COMM_WORLD, MPI_STATUS_IGNORE));

        localsum *= prevMultiplier;
        multiplier *= prevMultiplier;
    }

    if (rank < size - 1) {
        // Send multiplier to next rank
        PetscCall(MPI_Send(&multiplier, 1, MPIU_REAL, rank + 1, 0,
            PETSC_COMM_WORLD));
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
