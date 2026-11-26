#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "spm.h"
#include "fespm.h"

ErrorCode FESLSPMComputeTimestep_CFL(FESLSPM s, double uh[], double *dt);

/* FEM diffusion is explicit in time. The mass-matrix system is solved
 * at each step, and stability is controlled by FESPMGetStableTimestep(). */

ErrorCode FESLSPMComputeTimestep_const(FESLSPM s, double *dt)
{
  ErrorCode ierr;
  *dt = s->dt_const;
  return (0);
}

ErrorCode FESLSPMComputeTimestep(FESLSPM s, double h[], double *dt)
{
  ErrorCode ierr;

  switch (s->dtmode)
  {
  case DTM_CONST:
    ierr = FESLSPMComputeTimestep_const(s, dt);
    break;
  case DTM_EXPLICIT:
  {
    double dt_adv;
    double dt_diff = 1.0e32;

    /* Horizontal advection is constrained by its CFL condition. */
    ierr = FESLSPMComputeTimestep_CFL(s, s->u_h, &dt_adv);
    if (ierr != 0)
      return ierr;

    /* FEM diffusion is explicit in time, even though each step requires
       solving a mass-matrix linear system. Enforce its spectral stability
       limit as well as the advection CFL limit. */
    if (s->fespm != NULL)
    {
      ierr = FESPMGetStableTimestep(s->fespm, &dt_diff);
      if (ierr != 0)
      {
        ERRQ("Failed to compute FEM diffusion stability timestep");
        return ierr;
      }
    }

    *dt = (dt_adv < dt_diff) ? dt_adv : dt_diff;
    break;
  }
  default:
    ERRQ("Unknown timestep mode detected");
    return (1);
  }

  if (*dt < 0.0 || isnan(*dt) || isinf(*dt))
    ERRQ("Invalid timestep detected");
  /* For explicit stepping, dt_min must never override either the
     advection CFL bound or the FEM diffusion stability bound. */
  if (*dt < s->dt_min && s->dtmode == DTM_CONST)
    *dt = s->dt_min;
  if (*dt > s->dt_max)
    *dt = s->dt_max;
  if (*dt > s->time_end)
    ERRQ("Timestep is larger than final time");
  return (0);
}

ErrorCode FESLSPMUpdate_EnforceBoundaryConditions(FESLSPM s, double field[])
{
  int ni, nj, nidx_ij, nidx_neighbour;

  /* update interior first to avoid order dependence */
  if (s->bctype[DIR_JMAX] == 'N')
  { /* north */
    nj = s->ny - 1;
    for (ni = 1; ni < s->nx - 1; ni++)
    {
      nidx_ij = NodeIdx(s, ni, nj);
      nidx_neighbour = NodeIdx(s, ni, nj - 1);
      field[nidx_ij] = field[nidx_neighbour];
    }
  }
  if (s->bctype[DIR_JMIN] == 'N')
  { /* south */
    nj = 0;
    for (ni = 1; ni < s->nx - 1; ni++)
    {
      nidx_ij = NodeIdx(s, ni, nj);
      nidx_neighbour = NodeIdx(s, ni, nj + 1);
      field[nidx_ij] = field[nidx_neighbour];
    }
  }
  if (s->bctype[DIR_IMAX] == 'N')
  { /* east */
    ni = s->nx - 1;
    for (nj = 1; nj < s->ny - 1; nj++)
    {
      nidx_ij = NodeIdx(s, ni, nj);
      nidx_neighbour = NodeIdx(s, ni - 1, nj);
      field[nidx_ij] = field[nidx_neighbour];
    }
  }
  if (s->bctype[DIR_IMIN] == 'N')
  { /* west */
    ni = 0;
    for (nj = 1; nj < s->ny - 1; nj++)
    {
      nidx_ij = NodeIdx(s, ni, nj);
      nidx_neighbour = NodeIdx(s, ni + 1, nj);
      field[nidx_ij] = field[nidx_neighbour];
    }
  }

  /* corner updates */
  if ((s->bctype[DIR_JMAX] == 'N') && (s->bctype[DIR_IMAX] == 'N'))
  { // NE
    ni = s->nx - 1;
    nj = s->ny - 1;
    nidx_ij = NodeIdx(s, ni, nj);
    nidx_neighbour = NodeIdx(s, ni - 1, nj - 1);
    field[nidx_ij] = field[nidx_neighbour];
  }

  if ((s->bctype[DIR_JMIN] == 'N') && (s->bctype[DIR_IMAX] == 'N'))
  { // SE
    ni = s->nx - 1;
    nj = 0;
    nidx_ij = NodeIdx(s, ni, nj);
    nidx_neighbour = NodeIdx(s, ni - 1, nj + 1);
    field[nidx_ij] = field[nidx_neighbour];
  }

  if ((s->bctype[DIR_JMIN] == 'N') && (s->bctype[DIR_IMIN] == 'N'))
  { // SW
    ni = 0;
    nj = 0;
    nidx_ij = NodeIdx(s, ni, nj);
    nidx_neighbour = NodeIdx(s, ni + 1, nj + 1);
    field[nidx_ij] = field[nidx_neighbour];
  }

  if ((s->bctype[DIR_IMIN] == 'N') && (s->bctype[DIR_JMAX] == 'N'))
  { // NW
    ni = 0;
    nj = s->ny - 1;
    nidx_ij = NodeIdx(s, ni, nj);
    nidx_neighbour = NodeIdx(s, ni + 1, nj - 1);
    field[nidx_ij] = field[nidx_neighbour];
  }

  return (0);
}

ErrorCode FESLSPMView(FESLSPM s)
{
  printf("FESLSPM:\n");
  printf("  domain [%1.4e,%1.4e] x [%1.4e,%1.4e] \n", s->xrange[0], s->xrange[1], s->yrange[0], s->yrange[1]);
  printf("  nx %4d\n", s->nx);
  printf("  ny %4d\n", s->ny);
  printf("  time_end %1.4e\n", s->time_end);
  printf("  maxsteps %d\n", s->nt_max);
  printf("  output every %d\n", s->nt_io);

  printf("  timestep mode %d\n", (int)s->dtmode);
  switch (s->dtmode)
  {
  case DTM_CONST:
    printf("    dt_const %1.4e\n", s->dt_const);
    break;
  case DTM_EXPLICIT:
    printf("    dt_theta %1.4e\n", s->dt_theta);
    printf("    dt(%D) %1.4e\n", s->nt, s->dt);
    break;
  }

  printf("  diffusion mode %d\n", (int)s->mode);
  switch (s->mode)
  {
  case DM_CONST:
    printf("    k0 %1.4e\n", s->k0);
    break;
  case DM_NONLINEAR_CRITICALSLOPE:
    printf("    k0 %1.4e , m %1.4e Sc %1.4e\n", s->k0, s->d_exp, s->Sc);
    break;
  case DM_NONLINEAR_POWERLAW:
    printf("    k0 %1.4e , m %1.4e\n", s->k0, s->d_exp);
    break;
  case DM_UNDEFINED:
    break;
  }

  return (0);
}

ErrorCode FESLSPMStepInitialize(FESLSPM s)
{
  ErrorCode ierr;

  s->time = 0.0;
  s->nt = 0;

  ierr = FESLSPMUpdate_EnforceBoundaryConditions(s, s->h_ic);
  memcpy(s->h, s->h_ic, sizeof(double) * s->nnodes);
  memcpy(s->h_c, s->h_ic, sizeof(double) * s->nnodes); /* I set h_c so that the difference in height on the IC is zero */

  ierr = FESLSPMCheckpoint(s);

  return (0);
}

ErrorCode FESLSPMStepAdvDiff(FESLSPM s)
{
  ErrorCode ierr;
  int tk;

  printf("____ FESLSPM AdvDiff integrator ____\n");
  printf("  current time/step %1.4e/%d : target time %1.4e : max steps %d\n", s->time, s->nt, s->time_end, s->nt_max);

  if (s->time >= s->time_end)
  {
    printf("  FESLSPM: No updates performed as already at target end time\n");
    return (0);
  }

  for (tk = 0; tk < s->nt_max; tk++)
  {
    /* FEM assembles the element diffusivity from the current h* field. */
    ierr = FESLSPMComputeTimestep(s, s->h, &s->dt);

    if (s->time + s->dt > s->time_end)
    {
      s->dt = s->time_end - s->time;
    }

    /*
     * --------------------------------------------------------
     * 1. Existing Semi-Lagrangian advection
     * --------------------------------------------------------
     */

    ierr =
        FESLSPMUpdate_SemiLagrangian(
            s,
            s->dt,
            s->h);

    /*
     * h_c = h*
     *
     * h* is the result of the
     * Semi-Lagrangian advection step.
     */
    memcpy(
        s->h_c,
        s->h,
        sizeof(double) * s->nnodes);

    /*
     * --------------------------------------------------------
     * 2. FEM diffusion + vertical velocity
     * --------------------------------------------------------
     *
     * FESPM solves:
     *
     *     (h^(n+1) - h*) / dt
     *
     *        = u_v
     *          + div(kappa grad h*)
     *
     * --------------------------------------------------------
     */

    ierr =
        FESPMDiffusionStep(
            s->fespm,
            s,
            s->h_c);

    s->time += s->dt;
    s->nt++;
    printf("  step %4d: time %1.4e, dt %1.4e\n", s->nt, s->time, s->dt);

    if (s->nt % s->nt_io == 0)
    {
      ierr = FESLSPMFieldView(s, FVT_VTRB);
      ierr = FESLSPMCheckpoint(s);
    }
    if (s->time >= s->time_end)
    {
      printf("  FESLSPM: Reached target end time\n");
      break;
    }
  }
  if (tk == (s->nt_max))
  {
    printf("  FESLSPM: Reached max number of time steps permitted\n");
  }

  return (0);
}

/*
 Accepts user defined call back function: user_update_function() (may be NULL)
 The user function is called at each time step.
 The call back enables user to change velocity, boundary conditions as a function of time

 The prototype for the call back function is
   ErrorCode user_update_function(FESLSPM s,double time,int step,double dt,void *ctx)

 The pointer ctx is for any user-dependent data which is required by user_update_function().
*/

