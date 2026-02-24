//STARTWHOLE
static char help[] = "Solve a tridiagonal system of arbitrary size.\n"
"Option prefix = tri_.\n";

#include <petsc.h>
#include "petscviewerhdf5.h" 

int main(int argc,char **args) {
    Vec         u, f, uexact;
    Mat         A;
    KSP         ksp;
    PetscInt    m = 4, i, Istart, Iend, j[3];
    PetscReal   v[3], xval, errnorm;

    PetscReal   gamma = 0.0, c = 0.0, h, xi, rhs, uexactnorm, relerr;
    PetscInt    k = 1;
    PetscInt    rowsBC[2];

    PetscViewer  viewer;

    PetscCall(PetscInitialize(&argc,&args,NULL,help));

    PetscOptionsBegin(PETSC_COMM_WORLD,"bvp_","options for bvp",NULL);
    PetscCall(PetscOptionsInt("-m","dimension of linear system","bvp.c",m,&m,NULL));
    // BVP parameters
    PetscCall(PetscOptionsReal("-gamma","BVP gamma value","bvp.c",gamma,&gamma,NULL));
    PetscCall(PetscOptionsInt("-k","BVP k value","bvp.c",k,&k,NULL));
    PetscCall(PetscOptionsReal("-c","BVP c value","bvp.c",c,&c,NULL));
    PetscOptionsEnd();

    PetscCall(VecCreate(PETSC_COMM_WORLD,&u));
    PetscCall(VecSetSizes(u,PETSC_DECIDE,m));
    PetscCall(VecSetFromOptions(u));
    PetscCall(VecDuplicate(u,&f));
    PetscCall(VecDuplicate(u,&uexact));

    PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
    PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m));
    PetscCall(MatSetOptionsPrefix(A,"a_"));
    PetscCall(MatSetFromOptions(A));
    PetscCall(MatSetUp(A));
    PetscCall(MatGetOwnershipRange(A,&Istart,&Iend));

    h = 1.0 / (PetscReal)(m - 1);

    for (i=Istart; i<Iend; i++) {

        if (i == 0 || i == m-1) {
            // Fill diagonal. BCs enforced later
            v[0] = 1.0;
            j[0] = i;
            PetscCall(MatSetValues(A,1,&i,1,j,v,INSERT_VALUES));
        } else {
            v[0] = -1.0/(h*h);
            v[1] =  2.0/(h*h) + gamma;
            v[2] = -1.0/(h*h);
            j[0] = i-1; j[1] = i; j[2] = i+1;
            PetscCall(MatSetValues(A,1,&i,3,j,v,INSERT_VALUES));
        }

        // manufactured solution and RHS f(x)
        xi   = (PetscReal)i * h;
        xval = PetscSinReal((PetscReal)k * PETSC_PI * xi) + c * PetscPowReal(xi - 0.5, 3.0);

        rhs  = ((PetscReal)(k*k) * PETSC_PI * PETSC_PI + gamma) * PetscSinReal((PetscReal)k * PETSC_PI * xi)
             - 6.0 * c * (xi - 0.5)
             + gamma * c * PetscPowReal(xi - 0.5, 3.0);

        PetscCall(VecSetValues(uexact,1,&i,&xval,INSERT_VALUES));
        PetscCall(VecSetValues(f,     1,&i,&rhs, INSERT_VALUES));
    }

    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
    PetscCall(VecAssemblyBegin(uexact));
    PetscCall(VecAssemblyEnd(uexact));
    PetscCall(VecAssemblyBegin(f));
    PetscCall(VecAssemblyEnd(f));

    // Enforce Dirichlet BC
    rowsBC[0] = 0;
    rowsBC[1] = m-1;
    PetscCall(MatZeroRowsColumns(A,2,rowsBC,1.0,uexact,f));

    PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
    PetscCall(KSPSetOperators(ksp,A,A));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPSolve(ksp,f,u));

    // Compute relative error
    PetscCall(VecNorm(uexact,NORM_2,&uexactnorm));
    PetscCall(VecAXPY(u,-1.0,uexact));
    PetscCall(VecNorm(u,NORM_2,&errnorm));
    relerr = (uexactnorm > 0.0) ? (errnorm/uexactnorm) : errnorm;

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
    "relative error for m = %d system is ||u-uexact||_2/||uexact||_2 = %.1e\n",m,(double)relerr));

    // output the solution, rhs, and exact solution to an HDF5 file
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, "bvp_solution.h5",
        FILE_MODE_WRITE, &viewer));
    PetscCall(PetscObjectSetName((PetscObject) uexact, "uexact"));
    PetscCall(PetscObjectSetName((PetscObject) f, "f"));
    PetscCall(PetscObjectSetName((PetscObject) u, "u"));
    PetscCall(VecView(f, viewer));
    PetscCall(VecView(u, viewer));
    PetscCall(VecView(uexact, viewer));
    PetscCall(PetscViewerDestroy(&viewer));

    PetscCall(KSPDestroy(&ksp));
    PetscCall(MatDestroy(&A));
    PetscCall(VecDestroy(&u));
    PetscCall(VecDestroy(&f));
    PetscCall(VecDestroy(&uexact));
    PetscCall(PetscFinalize());
    return 0;
}
//ENDWHOLE