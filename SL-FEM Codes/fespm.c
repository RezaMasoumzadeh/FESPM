#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fespm.h"


/*
 * ============================================================
 * FESPM
 * ============================================================
 *
 * FEM diffusion module for the existing FESLSPM structured grid.
 *
 *
 * Mesh:
 *
 * Existing FESLSPM quadrilateral cells are divided into two
 * P1 triangular finite elements.
 *
 *
 *       n3 -------- n2
 *       |         / |
 *       |       /   |
 *       |     /     |
 *       |   /       |
 *       | /         |
 *       n0 -------- n1
 *
 *
 * Triangles:
 *
 *      (n0,n1,n2)
 *      (n0,n2,n3)
 *
 *
 * PDE:
 *
 *      dh/dt = u_v + div(kappa grad h)
 *
 *
 * After the Semi-Lagrangian step:
 *
 *      h_old = h*
 *
 *
 * Time discretization:
 *
 *      (H^{n+1} - H*) / dt
 *
 *          = U_v + div(kappa grad H*)
 *
 *
 * Galerkin FEM:
 *
 *      M (H^{n+1} - H*) / dt
 *
 *          + K H*
 *
 *          = F
 *
 *
 * where
 *
 *      M_ij = integral(phi_i phi_j) dOmega
 *
 *      K_ij =
 *          integral(
 *              kappa grad(phi_i).grad(phi_j)
 *          ) dOmega
 *
 *
 * If U_v is represented by nodal values:
 *
 *      U_v^h = sum_j U_v,j phi_j
 *
 * then
 *
 *      F = M U_v
 *
 *
 * Therefore:
 *
 *      M H^{n+1}
 *
 *        = M H*
 *
 *          + dt ( M U_v - K H* )
 *
 *
 * IMPORTANT:
 *
 * 1. Consistent mass matrix is used.
 *
 * 2. No mass lumping is used.
 *
 * 3. A linear system is solved at every timestep.
 *
 * 4. Diffusion is still evaluated explicitly at H*.
 *
 * 5. Therefore this is NOT backward Euler.
 *
 *
 * Boundary condition:
 *
 *      homogeneous Neumann
 *
 * is treated naturally through the weak formulation.
 *
 * Current version supports:
 *
 *      bctype = NNNN
 *
 * and constant or nonlinear critical-slope diffusivity.
 *
 * ============================================================
 */


/*
 * Maximum number of nonzero entries per row.
 *
 * For the structured P1 triangular mesh used here,
 * 9 entries are sufficient.
 */
#define FESPM_NNZ 9


/*
 * CG solver tolerance.
 */
#define FESPM_CG_TOL 1.0e-10


/*
 * Maximum CG iterations:
 *
 *      10 * number of nodes
 */
#define FESPM_CG_MAXIT_FACTOR 10


/*
 * Power iteration parameters used by
 * FESPMGetStableTimestep().
 */
#define FESPM_POWER_MAXIT 60
#define FESPM_POWER_TOL 1.0e-8


/*
 * ============================================================
 * FESPM data structure
 * ============================================================
 */

struct _p_FESPM {

    /* Parent surface solver, used to refresh nonlinear K for stability checks. */
    FESLSPM spm;

    int nx;
    int ny;

    int nnodes;
    int nelements;

    double dx;
    double dy;


    /*
     * --------------------------------------------------------
     * Consistent FEM mass matrix:
     *
     *      M
     *
     * stored using a fixed 9-entry sparse row.
     * --------------------------------------------------------
     */

    int    *M_col;
    double *M_val;


    /*
     * --------------------------------------------------------
     * FEM stiffness matrix:
     *
     *      K
     *
     * also stored using a fixed 9-entry sparse row.
     * --------------------------------------------------------
     */

    int    *K_col;
    double *K_val;


    /*
     * --------------------------------------------------------
     * Temporary vectors.
     * --------------------------------------------------------
     */


    /*
     *      MH = M H_old
     */
    double *MH;


    /*
     *      MU = M U_v
     */
    double *MU;


    /*
     *      KH = K H_old
     */
    double *KH;


    /*
     * Right-hand side:
     *
     *      RHS =
     *
     *      M H_old
     *
     *      + dt (M U_v - K H_old)
     */
    double *rhs;


    /*
     * New height:
     */
    double *Hnew;


    /*
     * --------------------------------------------------------
     * Conjugate Gradient work vectors.
     * --------------------------------------------------------
     */

    double *cg_r;
    double *cg_p;
    double *cg_Ap;
};


/*
 * ============================================================
 * Add an entry to a fixed 9-entry sparse row.
 * ============================================================
 */

static ErrorCode FESPMAddSparseEntry(
        int    *cols,
        double *vals,
        int row,
        int col,
        double value)
{
    int slot;


    /*
     * Search for an existing column.
     */
    for (slot = 0;
         slot < FESPM_NNZ;
         slot++) {

        if (cols[
                FESPM_NNZ * row + slot
            ] == col) {

            vals[
                FESPM_NNZ * row + slot
            ] += value;

            return(0);
        }
    }


    /*
     * Search for an empty slot.
     *
     * -1 means unused.
     */
    for (slot = 0;
         slot < FESPM_NNZ;
         slot++) {

        if (cols[
                FESPM_NNZ * row + slot
            ] == -1) {

            cols[
                FESPM_NNZ * row + slot
            ] = col;

            vals[
                FESPM_NNZ * row + slot
            ] = value;

            return(0);
        }
    }


    /*
     * More than 9 entries were required.
     */
    printf(
        "FESPM ERROR: more than %d entries "
        "in row %d.\n",
        FESPM_NNZ,
        row
    );

    return(1);
}


/*
 * ============================================================
 * Compute P1 triangle geometry.
 * ============================================================
 *
 * For a triangle:
 *
 *      grad(phi_a)
 *
 *          =
 *
 *      [ b_a ]
 *      [ c_a ]
 *      -------
 *        2A
 *
 * ============================================================
 */

static double FESPMTriangleGeometry(
        FESLSPM spm,
        int nodes[3],
        double x[3],
        double y[3],
        double b[3],
        double c[3])
{
    double area;


    /*
     * Coordinates of node 0.
     */
    x[0] =
        spm->xrange[0]
        + (nodes[0] % spm->nx)
        * spm->dx;

    y[0] =
        spm->yrange[0]
        + (nodes[0] / spm->nx)
        * spm->dy;


    /*
     * Coordinates of node 1.
     */
    x[1] =
        spm->xrange[0]
        + (nodes[1] % spm->nx)
        * spm->dx;

    y[1] =
        spm->yrange[0]
        + (nodes[1] / spm->nx)
        * spm->dy;


    /*
     * Coordinates of node 2.
     */
    x[2] =
        spm->xrange[0]
        + (nodes[2] % spm->nx)
        * spm->dx;

    y[2] =
        spm->yrange[0]
        + (nodes[2] / spm->nx)
        * spm->dy;


    /*
     * Triangle area.
     */
    area =
        0.5 *
        fabs(
            (x[1] - x[0])
            * (y[2] - y[0])
            -
            (x[2] - x[0])
            * (y[1] - y[0])
        );


    /*
     * P1 basis-function coefficients.
     */

    b[0] =
        y[1] - y[2];

    b[1] =
        y[2] - y[0];

    b[2] =
        y[0] - y[1];


    c[0] =
        x[2] - x[1];

    c[1] =
        x[0] - x[2];

    c[2] =
        x[1] - x[0];


    return(area);
}


/*
 * ============================================================
 * Assemble the CONSISTENT local mass matrix.
 * ============================================================
 *
 *
 * For a P1 triangle:
 *
 *
 *             A
 *      M^e = ------
 *             12
 *
 *
 *             [ 2  1  1 ]
 *
 *          x  [ 1  2  1 ]
 *
 *             [ 1  1  2 ]
 *
 *
 * This is NOT a lumped mass matrix.
 *
 * ============================================================
 */

static ErrorCode FESPMAddTriangleMass(
        FESPM fespm,
        int nodes[3],
        double area)
{
    static const double local[3][3] = {

        {
            2.0,
            1.0,
            1.0
        },

        {
            1.0,
            2.0,
            1.0
        },

        {
            1.0,
            1.0,
            2.0
        }
    };


    int a;
    int d;

    double factor;


    factor =
        area / 12.0;


    for (a = 0;
         a < 3;
         a++) {

        int row;

        row =
            nodes[a];


        for (d = 0;
             d < 3;
             d++) {

            int col;

            double value;


            col =
                nodes[d];


            value =
                factor
                * local[a][d];


            if (FESPMAddSparseEntry(
                    fespm->M_col,
                    fespm->M_val,
                    row,
                    col,
                    value)) {

                return(1);
            }
        }
    }


    return(0);
}


/*
 * ============================================================
 * Assemble local P1 stiffness matrix.
 * ============================================================
 *
 *
 *      K^e_ad
 *
 *        =
 *
 *      kappa
 *
 *      (b_a b_d + c_a c_d)
 *      --------------------
 *             4 A
 *
 * ============================================================
 */

static ErrorCode FESPMAddTriangleStiffness(
        FESPM fespm,
        int nodes[3],
        double b[3],
        double c[3],
        double area,
        double kappa)
{
    int a;
    int d;

    double factor;


    factor =
        kappa
        / (4.0 * area);


    for (a = 0;
         a < 3;
         a++) {

        int row;


        row =
            nodes[a];


        for (d = 0;
             d < 3;
             d++) {

            int col;

            double value;


            col =
                nodes[d];


            value =
                factor *
                (
                    b[a] * b[d]
                    +
                    c[a] * c[d]
                );


            if (FESPMAddSparseEntry(
                    fespm->K_col,
                    fespm->K_val,
                    row,
                    col,
                    value)) {

                return(1);
            }
        }
    }


    return(0);
}


/*
 * ============================================================
 * Assemble M and K.
 * ============================================================
 *
 * This is called once during FESPMCreate().
 *
 * ============================================================
 */

/*
 * ============================================================
 * Compute elementwise critical-slope diffusivity for a P1 triangle.
 *
 * For P1 elements, grad(h) is constant on each triangle, so
 * kappa is evaluated once per triangle for this constitutive law.
 *
 * The clipping at the critical slope follows the original FESLSPM
 * implementation: it prevents division by zero and caps kappa.
 * ============================================================
 */
static ErrorCode FESPMComputeElementKappa(
        FESLSPM spm,
        const int nodes[3],
        const double b[3],
        const double c[3],
        double area,
        const double h[],
        double *kappa_out)
{
    int a;
    double dhdx = 0.0;
    double dhdy = 0.0;
    double slope;
    double ratio;
    double nonlinear_term;
    double denominator;

    if (!spm || !nodes || !b || !c || !h || !kappa_out || area <= 0.0) {
        printf("FESPM ERROR: invalid arguments in element diffusivity.\n");
        return 1;
    }

    if (spm->k0 <= 0.0 || spm->d_exp <= 0.0) {
        printf("FESPM ERROR: k0 and d_exp must be positive.\n");
        return 1;
    }

    for (a = 0; a < 3; ++a) {
        dhdx += h[nodes[a]] * b[a] / (2.0 * area);
        dhdy += h[nodes[a]] * c[a] / (2.0 * area);
    }

    /* P1 triangular elements have a constant gradient per element. */
    slope = sqrt(dhdx * dhdx + dhdy * dhdy);

    if (spm->mode == DM_NONLINEAR_POWERLAW) {
        /* sfcmod=5: kappa = k0 * |grad(h)|^d_exp. */
        *kappa_out = spm->k0 * pow(slope, spm->d_exp);
    } else if (spm->mode == DM_NONLINEAR_CRITICALSLOPE) {
        /* sfcmod=6: kappa = k0 / (1 - (|grad(h)|/Sc)^d_exp). */
        if (spm->Sc <= 0.0) {
            printf("FESPM ERROR: Sc must be positive for critical-slope mode.\n");
            return 1;
        }
        ratio = slope / spm->Sc;
        if (ratio >= 1.0) {
            /* Match the original FDSPM regularization near the critical slope. */
            nonlinear_term = 1.0 - 1.0e-5;
        } else {
            nonlinear_term = pow(ratio, spm->d_exp);
        }
        denominator = 1.0 - nonlinear_term;
        if (!isfinite(denominator) || denominator <= 0.0) {
            printf("FESPM ERROR: invalid critical-slope denominator.\n");
            return 1;
        }
        *kappa_out = spm->k0 / denominator;
    } else {
        printf("FESPM ERROR: element diffusivity requested for an unsupported mode.\n");
        return 1;
    }

    /* Zero kappa is valid for a perfectly flat power-law element. */
    if (!isfinite(*kappa_out) || *kappa_out < 0.0) {
        printf("FESPM ERROR: invalid element diffusivity.\n");
        return 1;
    }

    return 0;
}


/* Reset only the sparse stiffness matrix before reassembly. */
static void FESPMClearStiffness(FESPM fespm)
{
    int i;
    int n = FESPM_NNZ * fespm->nnodes;

    for (i = 0; i < n; ++i) {
        fespm->K_col[i] = -1;
        fespm->K_val[i] = 0.0;
    }
}


/*
 * ============================================================
 * Assemble M and/or K.
 *
 * Parameters:
 *   assemble_mass      : build the consistent mass matrix
 *   assemble_stiffness : build the stiffness matrix
 *   store_kappa        : if nonzero, write the per-cell average of
 *                        the FEM per-triangle kappa into
 *                        spm->kappa[cell_idx]. This is the diagnostic
 *                        value reported in the .vtk/.vtr output and
 *                        checkpoint. Set to 1 only when this call
 *                        corresponds to the actual diffusion solve.
 * ============================================================
 */
static ErrorCode FESPMAssemble(
        FESPM fespm,
        FESLSPM spm,
        const double h[],
        int assemble_mass,
        int assemble_stiffness,
        int store_kappa)
{
    int i;
    int j;

    double kappa;


    if (spm->mode != DM_CONST &&
        spm->mode != DM_NONLINEAR_POWERLAW &&
        spm->mode != DM_NONLINEAR_CRITICALSLOPE) {
        printf("FESPM ERROR: unsupported diffusivity mode.\n");
        return 1;
    }

    if (!h) {
        printf("FESPM ERROR: surface field is required for assembly.\n");
        return 1;
    }

    if (spm->mode == DM_CONST) {
        kappa = spm->k0;
        if (kappa <= 0.0 || !isfinite(kappa)) {
            printf("FESPM ERROR: constant kappa must be positive and finite.\n");
            return 1;
        }
    } else if (spm->k0 <= 0.0 || spm->d_exp <= 0.0 ||
               (spm->mode == DM_NONLINEAR_CRITICALSLOPE && spm->Sc <= 0.0)) {
        printf("FESPM ERROR: invalid nonlinear diffusivity parameters.\n");
        return 1;
    }

    if (assemble_stiffness) {
        FESPMClearStiffness(fespm);
    }


    /*
     * Loop over FESLSPM quadrilateral cells.
     */

    for (j = 0;
         j < spm->my;
         j++) {

        for (i = 0;
             i < spm->mx;
             i++) {

            int n0;
            int n1;
            int n2;
            int n3;

            double kappa_tri1 = 0.0;
            double kappa_tri2 = 0.0;


            /*
             * Four corners of the quadrilateral.
             */

            n0 =
                i
                + j * spm->nx;


            n1 =
                (i + 1)
                + j * spm->nx;


            n2 =
                (i + 1)
                + (j + 1)
                * spm->nx;


            n3 =
                i
                + (j + 1)
                * spm->nx;


            /*
             * ------------------------------------------------
             * Triangle 1:
             *
             *      n2
             *      /|
             *     / |
             *    /  |
             *   /   |
             *  /    |
             * n0 --- n1
             *
             * (n0,n1,n2)
             * ------------------------------------------------
             */

            {
                int nodes[3];

                double x[3];
                double y[3];

                double b[3];
                double c[3];

                double area;


                nodes[0] =
                    n0;

                nodes[1] =
                    n1;

                nodes[2] =
                    n2;


                area =
                    FESPMTriangleGeometry(
                        spm,
                        nodes,
                        x,
                        y,
                        b,
                        c
                    );


                if (area <= 0.0) {

                    printf(
                        "FESPM ERROR: "
                        "invalid triangle area.\n"
                    );

                    return(1);
                }


                /*
                 * Consistent mass.
                 */
                if (assemble_mass &&
                    FESPMAddTriangleMass(
                        fespm,
                        nodes,
                        area)) {

                    return(1);
                }


                /*
                 * Stiffness.
                 */
                if (assemble_stiffness) {
                    double element_kappa = kappa;

                    if ((spm->mode == DM_NONLINEAR_POWERLAW ||
                         spm->mode == DM_NONLINEAR_CRITICALSLOPE) &&
                        FESPMComputeElementKappa(
                            spm, nodes, b, c, area, h,
                            &element_kappa)) {
                        return 1;
                    }

                    kappa_tri1 = element_kappa;

                    if (FESPMAddTriangleStiffness(
                            fespm,
                            nodes,
                            b,
                            c,
                            area,
                            element_kappa)) {
                        return 1;
                    }
                }
            }


            /*
             * ------------------------------------------------
             * Triangle 2:
             *
             *      n3 --- n2
             *       \   /
             *        \ /
             *         /
             *        /
             *       /
             *      n0
             *
             * (n0,n2,n3)
             * ------------------------------------------------
             */

            {
                int nodes[3];

                double x[3];
                double y[3];

                double b[3];
                double c[3];

                double area;


                nodes[0] =
                    n0;

                nodes[1] =
                    n2;

                nodes[2] =
                    n3;


                area =
                    FESPMTriangleGeometry(
                        spm,
                        nodes,
                        x,
                        y,
                        b,
                        c
                    );


                if (area <= 0.0) {

                    printf(
                        "FESPM ERROR: "
                        "invalid triangle area.\n"
                    );

                    return(1);
                }


                /*
                 * Consistent mass.
                 */
                if (assemble_mass &&
                    FESPMAddTriangleMass(
                        fespm,
                        nodes,
                        area)) {

                    return(1);
                }


                /*
                 * Stiffness.
                 */
                if (assemble_stiffness) {
                    double element_kappa = kappa;

                    if ((spm->mode == DM_NONLINEAR_POWERLAW ||
                         spm->mode == DM_NONLINEAR_CRITICALSLOPE) &&
                        FESPMComputeElementKappa(
                            spm, nodes, b, c, area, h,
                            &element_kappa)) {
                        return 1;
                    }

                    kappa_tri2 = element_kappa;

                    if (FESPMAddTriangleStiffness(
                            fespm,
                            nodes,
                            b,
                            c,
                            area,
                            element_kappa)) {
                        return 1;
                    }
                }
            }


            /*
             * Store the FEM-based per-cell kappa diagnostic.
             *
             * This is the value reported in the .vtk / .vtr output
             * and the checkpoint file. It reflects the exact kappa
             * that built the local stiffness matrices above.
             *
             * For DM_CONST, kappa_tri1 == kappa_tri2 == k0.
             * For DM_NONLINEAR_*, kappa_tri1 and kappa_tri2 come
             * from the P1 gradient of the current h field, at the
             * same time level as the actual solve (h*).
             */
            if (assemble_stiffness &&
                store_kappa &&
                spm->kappa != NULL) {

                spm->kappa[
                    i + j * spm->mx
                ] = 0.5 * (kappa_tri1 + kappa_tri2);
            }
        }
    }


    return(0);
}


/*
 * ============================================================
 * Sparse matrix-vector multiplication:
 *
 *      y = A x
 *
 * ============================================================
 */

static void FESPMMatVec(
        int *cols,
        double *vals,
        int n,
        const double x[],
        double y[])
{
    int i;
    int slot;


    for (i = 0;
         i < n;
         i++) {

        double sum;


        sum =
            0.0;


        for (slot = 0;
             slot < FESPM_NNZ;
             slot++) {

            int col;


            col =
                cols[
                    FESPM_NNZ * i
                    + slot
                ];


            if (col >= 0) {

                sum +=
                    vals[
                        FESPM_NNZ * i
                        + slot
                    ]
                    * x[col];
            }
        }


        y[i] =
            sum;
    }
}


/*
 * ============================================================
 * Dot product.
 * ============================================================
 */

static double FESPMDot(
        int n,
        const double a[],
        const double b[])
{
    int i;

    double sum;


    sum =
        0.0;


    for (i = 0;
         i < n;
         i++) {

        sum +=
            a[i]
            * b[i];
    }


    return(sum);
}


/*
 * ============================================================
 * Euclidean norm.
 * ============================================================
 */

static double FESPMNorm2(
        int n,
        const double a[])
{
    double value;


    value =
        FESPMDot(
            n,
            a,
            a
        );


    return(
        sqrt(value)
    );
}


/*
 * ============================================================
 * Conjugate Gradient solver.
 * ============================================================
 *
 * Solve:
 *
 *      M x = rhs
 *
 *
 * M is the consistent mass matrix.
 *
 * Since M is symmetric positive definite,
 * Conjugate Gradient is appropriate.
 *
 * ============================================================
 */

static ErrorCode FESPMCGSolveMass(
        FESPM fespm,
        const double rhs[],
        double x[])
{
    int n;

    int maxit;

    int i;


    double rhs_norm;

    double rr;

    double rr_new;


    n =
        fespm->nnodes;


    maxit =
        FESPM_CG_MAXIT_FACTOR
        * n;


    if (maxit < 100) {

        maxit =
            100;
    }


    /*
     * Initial guess:
     *
     *      x = 0
     */

    for (i = 0;
         i < n;
         i++) {

        x[i] =
            0.0;


        fespm->cg_r[i] =
            rhs[i];


        fespm->cg_p[i] =
            fespm->cg_r[i];


        fespm->cg_Ap[i] =
            0.0;
    }


    rhs_norm =
        FESPMNorm2(
            n,
            rhs
        );


    if (rhs_norm == 0.0) {

        return(0);
    }


    rr =
        FESPMDot(
            n,
            fespm->cg_r,
            fespm->cg_r
        );


    /*
     * CG iterations.
     */

    for (i = 0;
         i < maxit;
         i++) {

        double pAp;

        double alpha;

        double residual_norm;


        /*
         * Ap = M p
         */

        FESPMMatVec(
            fespm->M_col,
            fespm->M_val,
            n,
            fespm->cg_p,
            fespm->cg_Ap
        );


        /*
         * p^T A p
         */

        pAp =
            FESPMDot(
                n,
                fespm->cg_p,
                fespm->cg_Ap
            );


        if (pAp <= 0.0) {

            printf(
                "FESPM ERROR: CG breakdown: "
                "p^T M p <= 0.\n"
            );

            return(1);
        }


        /*
         * alpha
         */

        alpha =
            rr / pAp;


        /*
         * x = x + alpha p
         *
         * r = r - alpha Ap
         */

        for (int j = 0;
             j < n;
             j++) {

            x[j] +=
                alpha
                * fespm->cg_p[j];


            fespm->cg_r[j] -=
                alpha
                * fespm->cg_Ap[j];
        }


        /*
         * Check convergence.
         */

        residual_norm =
            FESPMNorm2(
                n,
                fespm->cg_r
            );


        if (residual_norm
            <=
            FESPM_CG_TOL
            * rhs_norm) {

            return(0);
        }


        /*
         * rr_new = r^T r
         */

        rr_new =
            FESPMDot(
                n,
                fespm->cg_r,
                fespm->cg_r
            );


        if (rr == 0.0) {

            return(0);
        }


        /*
         * beta
         */

        {
            double beta;


            beta =
                rr_new / rr;


            /*
             * p = r + beta p
             */

            for (int j = 0;
                 j < n;
                 j++) {

                fespm->cg_p[j] =
                    fespm->cg_r[j]
                    +
                    beta
                    * fespm->cg_p[j];
            }
        }


        rr =
            rr_new;
    }


    printf(
        "FESPM WARNING: CG did not converge "
        "within %d iterations.\n",
        maxit
    );


    return(1);
}


/*
 * ============================================================
 * FESPMCreate
 * ============================================================
 *
 * IMPORTANT:
 *
 * Call:
 *
 *      FESLSPMSetGeometry(...)
 *
 *      FESLSPMSetResolution(...)
 *
 *      FESLSPMSetConstantDiffusivity(...) OR
 *      FESLSPMSetNonlinearDiffusivityCriticalSlope(...)
 *
 *      FESPMCreate(...)
 *
 * in this order.
 *
 * ============================================================
 */

ErrorCode FESPMCreate(
        FESLSPM spm,
        FESPM *fespm_out)
{
    FESPM fespm;

    int n;


    if (!spm ||
        !fespm_out) {

        return(1);
    }


    *fespm_out =
        NULL;


    /* Constant and nonlinear critical-slope diffusivity are supported. */
    if (spm->mode != DM_CONST &&
        spm->mode != DM_NONLINEAR_POWERLAW &&
        spm->mode != DM_NONLINEAR_CRITICALSLOPE) {
        printf("FESPM ERROR: set a supported diffusivity mode before FESPMCreate().\n");
        return 1;
    }

    if ((spm->mode == DM_NONLINEAR_POWERLAW &&
         (spm->k0 <= 0.0 || spm->d_exp <= 0.0)) ||
        (spm->mode == DM_NONLINEAR_CRITICALSLOPE &&
         (spm->k0 <= 0.0 || spm->Sc <= 0.0 || spm->d_exp <= 0.0))) {
        printf("FESPM ERROR: invalid nonlinear diffusivity parameters.\n");
        return 1;
    }


    /*
     * Allocate FESPM structure.
     */

    fespm =
        (FESPM)
        malloc(
            sizeof(struct _p_FESPM)
        );


    if (!fespm) {

        printf(
            "FESPM ERROR: "
            "memory allocation failed.\n"
        );

        return(1);
    }


    memset(
        fespm,
        0,
        sizeof(struct _p_FESPM)
    );


    /*
     * Copy mesh information.
     */

    fespm->spm = spm;

    fespm->nx =
        spm->nx;

    fespm->ny =
        spm->ny;

    fespm->nnodes =
        spm->nnodes;

    fespm->nelements =
        2 * spm->ncells;

    fespm->dx =
        spm->dx;

    fespm->dy =
        spm->dy;


    /*
     * --------------------------------------------------------
     * Allocate consistent mass matrix.
     * --------------------------------------------------------
     */

    fespm->M_col =
        (int *)
        malloc(
            sizeof(int)
            * FESPM_NNZ
            * fespm->nnodes
        );


    fespm->M_val =
        (double *)
        malloc(
            sizeof(double)
            * FESPM_NNZ
            * fespm->nnodes
        );


    /*
     * --------------------------------------------------------
     * Allocate stiffness matrix.
     * --------------------------------------------------------
     */

    fespm->K_col =
        (int *)
        malloc(
            sizeof(int)
            * FESPM_NNZ
            * fespm->nnodes
        );


    fespm->K_val =
        (double *)
        malloc(
            sizeof(double)
            * FESPM_NNZ
            * fespm->nnodes
        );


    /*
     * --------------------------------------------------------
     * Allocate work vectors.
     * --------------------------------------------------------
     */

    fespm->MH =
        (double *)
        malloc(
            sizeof(double)
            * fespm->nnodes
        );


    fespm->MU =
        (double *)
        malloc(
            sizeof(double)
            * fespm->nnodes
        );


    fespm->KH =
        (double *)
        malloc(
            sizeof(double)
            * fespm->nnodes
        );


    fespm->rhs =
        (double *)
        malloc(
            sizeof(double)
            * fespm->nnodes
        );


    fespm->Hnew =
        (double *)
        malloc(
            sizeof(double)
            * fespm->nnodes
        );


    /*
     * --------------------------------------------------------
     * Allocate CG vectors.
     * --------------------------------------------------------
     */

    fespm->cg_r =
        (double *)
        malloc(
            sizeof(double)
            * fespm->nnodes
        );


    fespm->cg_p =
        (double *)
        malloc(
            sizeof(double)
            * fespm->nnodes
        );


    fespm->cg_Ap =
        (double *)
        malloc(
            sizeof(double)
            * fespm->nnodes
        );


    /*
     * Check allocations.
     */

    if (!fespm->M_col ||
        !fespm->M_val ||
        !fespm->K_col ||
        !fespm->K_val ||
        !fespm->MH ||
        !fespm->MU ||
        !fespm->KH ||
        !fespm->rhs ||
        !fespm->Hnew ||
        !fespm->cg_r ||
        !fespm->cg_p ||
        !fespm->cg_Ap) {

        printf(
            "FESPM ERROR: "
            "memory allocation failed.\n"
        );

        FESPMDestroy(
            &fespm
        );

        return(1);
    }


    /*
     * Initialize matrix values.
     */

    memset(
        fespm->M_val,
        0,
        sizeof(double)
        * FESPM_NNZ
        * fespm->nnodes
    );


    memset(
        fespm->K_val,
        0,
        sizeof(double)
        * FESPM_NNZ
        * fespm->nnodes
    );


    /*
     * -1 means unused sparse slot.
     */

    for (n = 0;
         n < FESPM_NNZ
             * fespm->nnodes;
         n++) {

        fespm->M_col[n] =
            -1;

        fespm->K_col[n] =
            -1;
    }


    /*
     * Assemble M and K once.
     *
     * store_kappa = 1 here: this is the initial diffusion state,
     * so spm->kappa[] should reflect the FEM kappa from h_init.
     */
    if (FESPMAssemble(
            fespm,
            spm,
            spm->h,
            1,  /* assemble mass matrix once */
            1,  /* assemble initial stiffness matrix */
            1   /* store FEM kappa into spm->kappa[] */
        )) {

        FESPMDestroy(
            &fespm
        );

        return(1);
    }


    *fespm_out =
        fespm;


    printf(
        "[FESPM] FEM diffusion module created.\n"
    );


    printf(
        "[FESPM] nodes = %d\n",
        fespm->nnodes
    );


    printf(
        "[FESPM] triangular elements = %d\n",
        fespm->nelements
    );


    printf(
        "[FESPM] consistent mass matrix M assembled.\n"
    );


    printf(
        "[FESPM] stiffness matrix K assembled.\n"
    );


    printf(
        "[FESPM] CG linear solve will be performed "
        "at every timestep.\n"
    );


    return(0);
}


/*
 * ============================================================
 * FESPMDestroy
 * ============================================================
 */

ErrorCode FESPMDestroy(
        FESPM *fespm)
{
    FESPM f;


    if (!fespm) {

        return(0);
    }


    f =
        *fespm;


    if (!f) {

        return(0);
    }


    /*
     * Matrices.
     */

    free(
        f->M_col
    );

    free(
        f->M_val
    );


    free(
        f->K_col
    );

    free(
        f->K_val
    );


    /*
     * Work vectors.
     */

    free(
        f->MH
    );

    free(
        f->MU
    );

    free(
        f->KH
    );

    free(
        f->rhs
    );

    free(
        f->Hnew
    );


    /*
     * CG vectors.
     */

    free(
        f->cg_r
    );

    free(
        f->cg_p
    );

    free(
        f->cg_Ap
    );


    free(
        f
    );


    *fespm =
        NULL;


    return(0);
}


/*
 * ============================================================
 * FESPMDiffusionStep
 * ============================================================
 *
 *
 * After Semi-Lagrangian advection:
 *
 *      h_old = h*
 *
 *
 * We solve:
 *
 *
 *      M Hnew
 *
 *        =
 *
 *      M H_old
 *
 *      + dt
 *
 *        ( M U_v
 *
 *          - K H_old )
 *
 *
 * ============================================================
 */

ErrorCode FESPMDiffusionStep(
        FESPM fespm,
        FESLSPM spm,
        double h_old[])
{
    int i;


    if (!fespm ||
        !spm ||
        !h_old) {

        return(1);
    }


    if (spm->mode != DM_CONST &&
        spm->mode != DM_NONLINEAR_POWERLAW &&
        spm->mode != DM_NONLINEAR_CRITICALSLOPE) {
        printf("FESPM ERROR: unsupported diffusivity mode.\n");
        return 1;
    }

    /*
     * Explicit diffusion: evaluate kappa from the known surface h_old.
     * Rebuild K only for the nonlinear mode; M remains unchanged.
     *
     * store_kappa = 1 here: this is the actual diffusion solve, so
     * spm->kappa[] must reflect the FEM kappa from h_old (= h*).
     */
    if (spm->mode == DM_NONLINEAR_POWERLAW ||
        spm->mode == DM_NONLINEAR_CRITICALSLOPE) {
        if (FESPMAssemble(
                fespm,
                spm,
                h_old,
                0,  /* do not reassemble mass matrix */
                1,  /* rebuild nonlinear stiffness matrix */
                1   /* store FEM kappa into spm->kappa[] */
            )) {
            printf("FESPM ERROR: nonlinear stiffness assembly failed.\n");
            return 1;
        }
    }


    /*
     * --------------------------------------------------------
     * Boundary conditions.
     *
     * Current implementation assumes:
     *
     *      NNNN
     *
     * homogeneous Neumann.
     *
     * These are natural BCs in the FEM weak form.
     * --------------------------------------------------------
     */

    for (i = 0;
         i < 4;
         i++) {

        if (spm->bctype[i] != 'N' &&
            spm->bctype[i] != 'n') {

            printf(
                "FESPM ERROR: current version "
                "supports only homogeneous "
                "Neumann boundary conditions (NNNN).\n"
            );

            return(1);
        }
    }


    /*
     * --------------------------------------------------------
     * Compute:
     *
     *      MH = M H_old
     * --------------------------------------------------------
     */

    FESPMMatVec(
        fespm->M_col,
        fespm->M_val,
        fespm->nnodes,
        h_old,
        fespm->MH
    );


    /*
     * --------------------------------------------------------
     * Compute:
     *
     *      MU = M U_v
     * --------------------------------------------------------
     */

    FESPMMatVec(
        fespm->M_col,
        fespm->M_val,
        fespm->nnodes,
        spm->u_v,
        fespm->MU
    );


    /*
     * --------------------------------------------------------
     * Compute:
     *
     *      KH = K H_old
     * --------------------------------------------------------
     */

    FESPMMatVec(
        fespm->K_col,
        fespm->K_val,
        fespm->nnodes,
        h_old,
        fespm->KH
    );


    /*
     * --------------------------------------------------------
     * Construct RHS:
     *
     *
     *      RHS
     *
     *        =
     *
     *      M H_old
     *
     *      + dt
     *
     *      ( M U_v - K H_old )
     *
     * --------------------------------------------------------
     */

    for (i = 0;
         i < fespm->nnodes;
         i++) {

        fespm->rhs[i] =

            fespm->MH[i]

            +

            spm->dt
            *
            (
                fespm->MU[i]
                -
                fespm->KH[i]
            );
    }


    /*
     * --------------------------------------------------------
     * Solve:
     *
     *      M Hnew = RHS
     *
     * --------------------------------------------------------
     */

    if (FESPMCGSolveMass(
            fespm,
            fespm->rhs,
            fespm->Hnew)) {

        printf(
            "FESPM ERROR: "
            "mass-matrix linear solve failed.\n"
        );

        return(1);
    }


    /*
     * --------------------------------------------------------
     * Copy FEM solution back into FESLSPM.
     * --------------------------------------------------------
     */

    memcpy(
        spm->h,
        fespm->Hnew,
        sizeof(double)
        * fespm->nnodes
    );


    /*
     * IMPORTANT:
     *
     * Do NOT call:
     *
     *      FESLSPMUpdate_EnforceBoundaryConditions()
     *
     * here.
     *
     * Homogeneous Neumann BCs are natural BCs in
     * the FEM weak formulation.
     */


    return(0);
}


/*
 * ============================================================
 * FESPMGetStableTimestep
 * ============================================================
 *
 * Because the present time discretization is explicit with
 * respect to the diffusion operator:
 *
 *
 *      M dH/dt + K H = F
 *
 *
 * the stability condition is approximately:
 *
 *
 *      dt <= 2 / lambda_max
 *
 *
 * where:
 *
 *
 *      K q = lambda M q
 *
 *
 * We estimate lambda_max using power iteration.
 *
 * ============================================================
 */

ErrorCode FESPMGetStableTimestep(
        FESPM fespm,
        double *dt)
{
    int n;

    int it;

    double *q;

    double *z;

    double *Kq;

    double lambda_old;

    double lambda;

    double qnorm;


    if (!fespm ||
        !dt) {

        return(1);
    }


    /* Refresh the nonlinear stiffness matrix from the current surface field. */
    if (fespm->spm &&
        (fespm->spm->mode == DM_NONLINEAR_POWERLAW ||
         fespm->spm->mode == DM_NONLINEAR_CRITICALSLOPE)) {
        /* store_kappa = 0: stability probing must not overwrite the
         * spm->kappa[] diagnostic that the diffusion solve wrote. */
        if (FESPMAssemble(fespm, fespm->spm, fespm->spm->h,
                          0, 1, 0)) {
            printf("FESPM ERROR: failed to refresh nonlinear K for stability estimate.\n");
            return 1;
        }
    }

    n =
        fespm->nnodes;


    q =
        (double *)
        malloc(
            sizeof(double)
            * n
        );


    z =
        (double *)
        malloc(
            sizeof(double)
            * n
        );


    Kq =
        (double *)
        malloc(
            sizeof(double)
            * n
        );


    if (!q ||
        !z ||
        !Kq) {

        free(q);
        free(z);
        free(Kq);

        return(1);
    }


    /*
     * Initial vector.
     *
     * Do not use a constant vector because K has
     * a zero eigenvalue for pure Neumann BCs.
     */

    for (int i = 0;
         i < n;
         i++) {

        q[i] =
            (i % 2 == 0)
            ? 1.0
            : -1.0;
    }


    qnorm =
        FESPMNorm2(
            n,
            q
        );


    for (int i = 0;
         i < n;
         i++) {

        q[i] /=
            qnorm;
    }


    lambda_old =
        0.0;


    lambda =
        0.0;


    /*
     * Power iteration.
     */

    for (it = 0;
         it < FESPM_POWER_MAXIT;
         it++) {

        double qMq;

        double qKq;


        /*
         * Kq = K q
         */

        FESPMMatVec(
            fespm->K_col,
            fespm->K_val,
            n,
            q,
            Kq
        );


        /*
         * z = M^{-1} Kq
         *
         * This requires a linear solve with M.
         */

        if (FESPMCGSolveMass(
                fespm,
                Kq,
                z)) {

            free(q);
            free(z);
            free(Kq);

            return(1);
        }


        /*
         * Mq.
         */

        FESPMMatVec(
            fespm->M_col,
            fespm->M_val,
            n,
            q,
            fespm->MH
        );


        /*
         * Rayleigh quotient:
         *
         *
         *      lambda
         *
         *       =
         *
         *      q^T K q
         *      ---------
         *      q^T M q
         *
         */

        qMq =
            FESPMDot(
                n,
                q,
                fespm->MH
            );


        qKq =
            FESPMDot(
                n,
                q,
                Kq
            );


        if (qMq <= 0.0) {

            free(q);
            free(z);
            free(Kq);

            return(1);
        }


        lambda =
            qKq
            / qMq;


        /*
         * Normalize z.
         */

        qnorm =
            FESPMNorm2(
                n,
                z
            );


        if (qnorm == 0.0) {

            break;
        }


        for (int i = 0;
             i < n;
             i++) {

            q[i] =
                z[i]
                / qnorm;
        }


        /*
         * Check convergence.
         */

        if (it > 0) {

            if (
                fabs(
                    lambda
                    -
                    lambda_old
                )
                <=
                FESPM_POWER_TOL
                *
                fabs(lambda)
            ) {

                break;
            }
        }


        lambda_old =
            lambda;
    }


    free(q);
    free(z);
    free(Kq);


    /*
     * If there is no positive eigenvalue,
     * return a very large timestep.
     */

    if (lambda <= 0.0) {

        *dt =
            1.0e32;

        return(0);
    }


    /*
     * Explicit Euler:
     *
     *      dt <= 2/lambda_max
     *
     *
     * Apply safety factor 0.9.
     */

    *dt =
        0.9
        * 2.0
        / lambda;


    return(0);
}