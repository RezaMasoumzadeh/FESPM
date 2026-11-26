#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "spm.h"
#include "fespm.h"

ErrorCode FESLSPMCreate(FESLSPM *s)
{
  ErrorCode ierr;
  FESLSPM spm;

  spm = malloc(sizeof(struct _p_FESLSPM));
  memset(spm, 0, sizeof(struct _p_FESLSPM));

  spm->dt_min = 1.0e-32;
  spm->dt_max = 1.0e32;
  spm->dt = -1.0;
  spm->dt_theta = 1.0;
  spm->mode = DM_CONST;
  spm->dtmode = DTM_EXPLICIT;

  spm->nt_io = 1;
  spm->nt_max = 1;
  spm->time_end = 1.0e32;

  spm->bctype[0] = 'D';
  spm->bctype[1] = 'D';
  spm->bctype[2] = 'D';
  spm->bctype[3] = 'D';
  spm->semilagrange_interpolation_type = SLINTERPOLATION_CSPLINE;

  *s = spm;
  return (0);
}

ErrorCode FESLSPMDestroy(FESLSPM *s)
{
  ErrorCode ierr;
  FESLSPM spm;

  if (!s)
    return (0);

  spm = *s;

  if (spm->h)
    free(spm->h);
  if (spm->h_c)
    free(spm->h_c);
  if (spm->h_ic)
    free(spm->h_ic);
  if (spm->u_v)
    free(spm->u_v);
  if (spm->u_h)
    free(spm->u_h);
  if (spm->k0_n)       /* ← h → k0_n */
    free(spm->k0_n);
  if (spm->d_exp_n)    /* ← h → d_exp_n */
    free(spm->d_exp_n);
  if (spm->Sc_n)       /* ← h → Sc_n */
    free(spm->Sc_n);
  if (spm->kappa)
    free(spm->kappa);

  if (spm->fespm)
  {
    FESPMDestroy(&spm->fespm);
  }
  free(spm);

  *s = NULL;

  return (0);
}

/* Mesh specification */
ErrorCode FESLSPMSetGeometry(FESLSPM s, double x0, double x1, double y0, double y1)
{
  ErrorCode ierr;

  s->xrange[0] = x0;
  s->xrange[1] = x1;
  s->yrange[0] = y0;
  s->yrange[1] = y1;

  return (0);
}

ErrorCode FESLSPMSetResolution(FESLSPM s, int nx, int ny)
{
  s->nx = nx;
  s->ny = ny;
  s->nnodes = nx * ny;
  s->mx = nx - 1;
  s->my = ny - 1;
  s->ncells = s->mx * s->my;

  s->dx = (s->xrange[1] - s->xrange[0]) / ((double)s->mx);
  s->dy = (s->yrange[1] - s->yrange[0]) / ((double)s->my);

  s->h = malloc(sizeof(double) * s->nnodes);
  memset(s->h, 0, sizeof(double) * s->nnodes);
  s->h_c = malloc(sizeof(double) * s->nnodes);
  memset(s->h_c, 0, sizeof(double) * s->nnodes);
  s->h_ic = malloc(sizeof(double) * s->nnodes);
  memset(s->h_ic, 0, sizeof(double) * s->nnodes);

  s->u_v = malloc(sizeof(double) * s->nnodes);
  memset(s->u_v, 0, sizeof(double) * s->nnodes);
  s->u_h = malloc(sizeof(double) * s->nnodes * 2);
  memset(s->u_h, 0, sizeof(double) * s->nnodes * 2);

  s->k0_n = malloc(sizeof(double) * s->nnodes);
  memset(s->k0_n, 0, sizeof(double) * s->nnodes);
  s->d_exp_n = malloc(sizeof(double) * s->nnodes);
  memset(s->d_exp_n, 0, sizeof(double) * s->nnodes);
  s->Sc_n = malloc(sizeof(double) * s->nnodes);
  memset(s->Sc_n, 0, sizeof(double) * s->nnodes);

  s->kappa =
      malloc(sizeof(double) * s->ncells);

  memset(
      s->kappa,
      0,
      sizeof(double) * s->ncells);

  return (0);
}

/* Timestepping specification */
ErrorCode FESLSPMSetTimestepBounds(FESLSPM spm, double dt_min, double dt_max)
{
  ErrorCode ierr;

  if (!FESLSPMCheckIsIgnore(dt_min))
  {
    spm->dt_min = dt_min;
  }
  if (!FESLSPMCheckIsIgnore(dt_max))
  {
    spm->dt_max = dt_max;
  }
  return (0);
}

ErrorCode FESLSPMSetTimestep(FESLSPM s, double dt)
{
  s->dt = dt;
  return (0);
}

ErrorCode FESLSPMSetTimestepTheta(FESLSPM s, double theta)
{
  s->dt_theta = theta;
  return (0);
}

ErrorCode FESLSPMSetTimestepIO(FESLSPM s, int nt)
{
  s->nt_io = nt;
  return (0);
}

ErrorCode FESLSPMSetMaxNumberSteps(FESLSPM s, int nt_max)
{
  s->nt_max = nt_max;
  return (0);
}

ErrorCode FESLSPMSetEndTime(FESLSPM s, double tend)
{
  s->time_end = tend;
  return (0);
}

ErrorCode FESLSPMSetTimestepMode(FESLSPM s, FESLSPMTimestepMode mode)
{
  s->dtmode = mode;
  return (0);
}

ErrorCode FESLSPMSetConstantTimestep(FESLSPM s, double dt)
{
  int ierr;
  ierr = FESLSPMSetTimestepMode(s, DTM_CONST);
  ierr = FESLSPMSetTimestep(s, dt);
  return (0);
}

ErrorCode FESLSPMSetExplicitTimestep(FESLSPM s, double theta)
{
  int ierr;
  ierr = FESLSPMSetTimestepMode(s, DTM_EXPLICIT);
  ierr = FESLSPMSetTimestepTheta(s, theta);
  return (0);
}

/* Diffusion specification */
ErrorCode FESLSPMSetDiffusionMode(FESLSPM s, FESLSPMDiffusionMode mode)
{
  s->mode = mode;
  return (0);
}

ErrorCode FESLSPMSetConstantDiffusivity(FESLSPM s, double k0)
{
  ErrorCode ierr;
  ierr = FESLSPMSetDiffusionMode(s, DM_CONST);
  s->k0 = k0;

  return (0);
}

/*
 slope angle = tan^{-1}(Sc)
*/
ErrorCode FESLSPMSetNonlinearDiffusivityCriticalSlope(FESLSPM s, double k0, double Sc, double d_exp)
{
  ErrorCode ierr;

  ierr = FESLSPMSetDiffusionMode(s, DM_NONLINEAR_CRITICALSLOPE);
  s->k0 = k0;
  s->d_exp = d_exp;
  s->Sc = Sc;

  return (0);
}

ErrorCode FESLSPMSetNonlinearDiffusivityPowerLaw(FESLSPM s, double k0, double d_exp)
{
  ErrorCode ierr;

  ierr = FESLSPMSetDiffusionMode(s, DM_NONLINEAR_POWERLAW);
  s->k0 = k0;
  s->d_exp = d_exp;

  return (0);
}

/* Initial condition helpers */
ErrorCode FESLSPMSetInitialHeight(FESLSPM s, double h0[])
{
  memcpy(s->h_ic, h0, sizeof(double) * s->nnodes);
  return (0);
}

ErrorCode FESLSPMSetVerticalVelocity(FESLSPM s, double v[])
{
  memcpy(s->u_v, v, sizeof(double) * s->nnodes);
  return (0);
}

ErrorCode FESLSPMSetHorizontalVelocity(FESLSPM s, double uw[])
{
  memcpy(s->u_h, uw, sizeof(double) * 2 * s->nnodes);
  return (0);
}

ErrorCode FESLSPMSetBoundaryConditions(FESLSPM s, const char bctype[])
{
  int k;

  memcpy(s->bctype, bctype, sizeof(char) * 4);
  for (k = 0; k < 4; k++)
  {
    char tc;

    tc = bctype[k];
    s->bctype[k] = toupper(tc);
    if (s->bctype[k] != 'D' && s->bctype[k] != 'N')
    {
      ERRQ("BCType must be N or D");
    }
  }
  return (0);
}

/* Mesh geometry helpers */
ErrorCode FESLSPMGetNodeCoordinate(FESLSPM s, int i, int j, double coor[])
{
  ErrorCode ierr;

  if (i < 0)
  {
    ERRQ("NodeCoordinate i < 0");
  }
  if (i >= s->nx)
  {
    ERRQ("NodeCoordinate i >= nx");
  }
  if (j < 0)
  {
    ERRQ("NodeCoordinate j < 0");
  }
  if (j >= s->ny)
  {
    ERRQ("NodeCoordinate j >= ny");
  }

  coor[0] = s->xrange[0] + i * s->dx;
  coor[1] = s->yrange[0] + j * s->dy;
  return (0);
}

ErrorCode FESLSPMGetCentroidCoordinate(FESLSPM s, int i, int j, double coor[])
{
  ErrorCode ierr;

  if (i < 0)
  {
    ERRQ("CellCoordinate i < 0");
  }
  if (i >= s->mx)
  {
    ERRQ("CellCoordinate i >= mx");
  }
  if (j < 0)
  {
    ERRQ("CellCoordinate j < 0");
  }
  if (j >= s->my)
  {
    ERRQ("CellCoordinate j >= my");
  }

  coor[0] = s->xrange[0] + 0.5 * s->dx + i * s->dx;
  coor[1] = s->yrange[0] + 0.5 * s->dy + j * s->dy;
  return (0);
}

ErrorCode FESLSPMFieldView_ascii(FESLSPM s)
{
  FILE *fp = NULL;
  char name[128];
  int nidx, cidx;
  double xn[2];

  sprintf(name, "spm-n-step%d.gp", s->nt);
  fp = fopen(name, "w");
  if (!fp)
    ERRQ2("Unable to open file", name);

  fprintf(fp, "# x (1) y (2) h (3) delta(h) (4) uz (5) ux (6) uy (7)\n");

  NodalIteratorBegin(s)
      nidx = NodeIdx(s, ni, nj);
  FESLSPMGetNodeCoordinate(s, ni, nj, xn);

  fprintf(fp, "%1.4e %1.4e ", xn[0], xn[1]);

  fprintf(fp, "%1.4e %1.4e ", s->h[nidx], s->h[nidx] - s->h_c[nidx]);

  fprintf(fp, "%1.4e %1.4e %1.4e\n", s->u_v[nidx], s->u_h[2 * nidx], s->u_h[2 * nidx + 1]);
  NodalIteratorEnd(s)
      fclose(fp);

  sprintf(name, "spm-c-step%d.gp", s->nt);
  fp = fopen(name, "w");
  if (!fp)
    ERRQ2("Unable to open file", name);
  CentroidIteratorBegin(s)
      cidx = CellIdx(s, mi, mj);
  FESLSPMGetCentroidCoordinate(s, mi, mj, xn);

  fprintf(fp, "%1.4e %1.4e ", xn[0], xn[1]);

  fprintf(fp, "%1.4e\n", s->kappa[cidx]);

  CentroidIteratorEnd(s)
      fclose(fp);

  return (0);
}

ErrorCode FESLSPMView_Binary(FESLSPM s, const char name[])
{
  ErrorCode ierr;
  FILE *fp;

  fp = fopen(name, "w");
  if (!fp)
    ERRQ2("Unable to open file", name);

  fwrite(s->xrange, sizeof(double), 2, fp);
  fwrite(s->yrange, sizeof(double), 2, fp);

  fwrite(&s->time_end, sizeof(double), 1, fp);
  fwrite(&s->dt, sizeof(double), 1, fp);
  fwrite(&s->dt_const, sizeof(double), 1, fp);
  fwrite(&s->time, sizeof(double), 1, fp);

  fwrite(&s->nx, sizeof(int), 1, fp);
  fwrite(&s->ny, sizeof(int), 1, fp);
  fwrite(&s->mx, sizeof(int), 1, fp);
  fwrite(&s->my, sizeof(int), 1, fp);
  fwrite(&s->nnodes, sizeof(int), 1, fp);
  fwrite(&s->ncells, sizeof(int), 1, fp);

  fwrite(&s->dx, sizeof(double), 1, fp);
  fwrite(&s->dy, sizeof(double), 1, fp);
  fwrite(&s->dt_theta, sizeof(double), 1, fp);

  fwrite(&s->nt, sizeof(int), 1, fp);
  fwrite(&s->nt_io, sizeof(int), 1, fp);
  fwrite(&s->nt_max, sizeof(int), 1, fp);

  fwrite(&s->k0, sizeof(double), 1, fp);
  fwrite(&s->d_exp, sizeof(double), 1, fp);
  fwrite(&s->Sc, sizeof(double), 1, fp);

  {
    int dummy;

    dummy = (int)s->mode;
    fwrite(&dummy, sizeof(int), 1, fp);
    dummy = (int)s->dtmode;
    fwrite(&dummy, sizeof(int), 1, fp);
  }

  fwrite(s->h_ic, sizeof(double), s->nnodes, fp);
  fwrite(s->h, sizeof(double), s->nnodes, fp);
  fwrite(s->h_c, sizeof(double), s->nnodes, fp);
  fwrite(s->u_v, sizeof(double), s->nnodes, fp);
  fwrite(s->u_h, sizeof(double), 2 * s->nnodes, fp);
  fwrite(s->kappa, sizeof(double), s->ncells, fp);
  fwrite(s->k0_n, sizeof(double), s->nnodes, fp);
  fwrite(s->d_exp_n, sizeof(double), s->nnodes, fp);
  fwrite(s->Sc_n, sizeof(double), s->nnodes, fp);
  fwrite(s->bctype, sizeof(char), 4, fp);

  fwrite(&s->dt_min, sizeof(double), 1, fp);
  fwrite(&s->dt_max, sizeof(double), 1, fp);

  fclose(fp);
  return (0);
}

ErrorCode FESLSPMFieldView_VTRB(FESLSPM s, const char filename[])
{
  FILE *fp = NULL;
  int i, offset, bytes;
  double coor;
  /* customization */
  int write_dh = 1;
  int write_vel = 1;
  int write_sedrate = 1;
  int write_kappa = 1;
  int write_k0 = 1;
  int write_dexp = 0;
  int write_sc = 0;

  fp = fopen(filename, "w");
  if (!fp)
    ERRQ2("Unable to open file", filename);

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  // fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");

  fprintf(fp, "  <RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, s->nx - 1, 0, s->ny - 1, 0, 0);
  fprintf(fp, "    <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, s->nx - 1, 0, s->ny - 1, 0, 0);

  offset = 0;

  fprintf(fp, "      <Coordinates>\n");
  fprintf(fp, "        <DataArray Name=\"cx\" type=\"Float64\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(double) * s->nx;
  fprintf(fp, "        <DataArray Name=\"cy\" type=\"Float64\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(double) * s->ny;
  fprintf(fp, "        <DataArray Name=\"cz\" type=\"Float64\" format=\"appended\" offset=\"%d\"/>\n", offset);
  offset += sizeof(int) + sizeof(double);
  fprintf(fp, "      </Coordinates>\n");

  fprintf(fp, "      <CellData>\n");
  if (write_kappa)
  {
    fprintf(fp, "        <DataArray Name=\"kappa_eff\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", offset);
    offset += sizeof(int) + sizeof(double) * s->ncells;
  }
  fprintf(fp, "      </CellData>\n");

  fprintf(fp, "      <PointData>\n");
  {
    fprintf(fp, "        <DataArray Name=\"h\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", offset);
    offset += sizeof(int) + sizeof(double) * s->nnodes;
  }
  if (write_dh)
  {
    fprintf(fp, "        <DataArray Name=\"dh\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", offset);
    offset += sizeof(int) + sizeof(double) * s->nnodes;
  }
  if (write_sedrate)
  {
    fprintf(fp, "        <DataArray Name=\"sed_rate\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", offset);
    offset += sizeof(int) + sizeof(double) * s->nnodes;
  }
  if (write_vel)
  {
    fprintf(fp, "        <DataArray Name=\"u_vert\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", offset);
    offset += sizeof(int) + sizeof(double) * s->nnodes;
    fprintf(fp, "        <DataArray Name=\"u_horiz\" type=\"Float64\" NumberOfComponents=\"2\" format=\"appended\" offset=\"%d\"/>\n", offset);
    offset += sizeof(int) + sizeof(double) * s->nnodes * 2;
  }
  if (write_k0)
  {
    fprintf(fp, "        <DataArray Name=\"k0\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", offset);
    offset += sizeof(int) + sizeof(double) * s->nnodes;
  }
  if (write_dexp)
  {
    if ((s->mode == DM_NONLINEAR_CRITICALSLOPE) || (s->mode == DM_NONLINEAR_POWERLAW))
    {
      fprintf(fp, "        <DataArray Name=\"dexp\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", offset);
      offset += sizeof(int) + sizeof(double) * s->nnodes;
    }
  }
  if (write_sc)
  {
    if (s->mode == DM_NONLINEAR_CRITICALSLOPE)
    {
      fprintf(fp, "        <DataArray Name=\"Sc\" type=\"Float64\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n", offset);
      offset += sizeof(int) + sizeof(double) * s->nnodes;
    }
  }

  fprintf(fp, "      </PointData>\n");

  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </RectilinearGrid>\n");

  fprintf(fp, "  <AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");

  /* <coordinates> x coordinates, y coordinates, dummy z coordinate */
  bytes = sizeof(double) * s->nx;
  fwrite(&bytes, sizeof(int), 1, fp);
  for (i = 0; i < s->nx; i++)
  {
    coor = s->xrange[0] + i * s->dx;
    fwrite(&coor, sizeof(double), 1, fp);
  }

  bytes = sizeof(double) * s->ny;
  fwrite(&bytes, sizeof(int), 1, fp);
  for (i = 0; i < s->ny; i++)
  {
    coor = s->yrange[0] + i * s->dy;
    fwrite(&coor, sizeof(double), 1, fp);
  }

  bytes = sizeof(double) * 1;
  fwrite(&bytes, sizeof(int), 1, fp);
  coor = 0.0;
  fwrite(&coor, sizeof(double), 1, fp);

  /* <CellData> kappa */
  if (write_kappa)
  {
    bytes = sizeof(double) * s->ncells;
    fwrite(&bytes, sizeof(int), 1, fp);
    fwrite(s->kappa, sizeof(double), s->ncells, fp);
  }

  /* <NodeData> h, delta_h, sed_rate, u_vertical, u_horiz */
  bytes = sizeof(double) * s->nnodes;
  fwrite(&bytes, sizeof(int), 1, fp);
  fwrite(s->h, sizeof(double), s->nnodes, fp);

  if (write_dh)
  {
    bytes = sizeof(double) * s->nnodes;
    fwrite(&bytes, sizeof(int), 1, fp);
    for (i = 0; i < s->nnodes; i++)
    {
      double delta_h;
      delta_h = s->h[i] - s->h_c[i];
      fwrite(&delta_h, sizeof(double), 1, fp);
    }
  }

  if (write_sedrate)
  {
    bytes = sizeof(double) * s->nnodes;
    fwrite(&bytes, sizeof(int), 1, fp);
    for (i = 0; i < s->nnodes; i++)
    {
      double sed;
      sed = (s->h[i] - s->h_c[i] - s->dt * s->u_v[i]) / s->dt;
      fwrite(&sed, sizeof(double), 1, fp);
    }
  }

  if (write_vel)
  {
    bytes = sizeof(double) * s->nnodes;
    fwrite(&bytes, sizeof(int), 1, fp);
    fwrite(s->u_v, sizeof(double), s->nnodes, fp);

    bytes = sizeof(double) * s->nnodes * 2;
    fwrite(&bytes, sizeof(int), 1, fp);
    fwrite(s->u_h, sizeof(double), s->nnodes * 2, fp);
  }

  if (write_k0)
  {
    bytes = sizeof(double) * s->nnodes;
    fwrite(&bytes, sizeof(int), 1, fp);
    fwrite(s->k0_n, sizeof(double), s->nnodes, fp);
  }

  if (write_dexp)
  {
    if ((s->mode == DM_NONLINEAR_CRITICALSLOPE) || (s->mode == DM_NONLINEAR_POWERLAW))
    {
      bytes = sizeof(double) * s->nnodes;
      fwrite(&bytes, sizeof(int), 1, fp);
      fwrite(s->d_exp_n, sizeof(double), s->nnodes, fp);
    }
  }

  if (write_sc)
  {
    if (s->mode == DM_NONLINEAR_CRITICALSLOPE)
    {
      bytes = sizeof(double) * s->nnodes;
      fwrite(&bytes, sizeof(int), 1, fp);
      fwrite(s->Sc_n, sizeof(double), s->nnodes, fp);
    }
  }

  fprintf(fp, "\n  </AppendedData>\n");

  fprintf(fp, "</VTKFile>\n");

  fclose(fp);

  return (0);
}

ErrorCode FESLSPMFieldView_VTR(FESLSPM s, const char filename[])
{
  FILE *fp = NULL;
  int i;
  double coor;

  fp = fopen(filename, "w");
  if (!fp)
    ERRQ2("Unable to open file", filename);

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  // fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");

  fprintf(fp, "  <RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, s->nx - 1, 0, s->ny - 1, 0, 0);
  fprintf(fp, "    <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, s->nx - 1, 0, s->ny - 1, 0, 0);

  fprintf(fp, "      <Coordinates>\n");
  fprintf(fp, "      <DataArray Name=\"cx\" type=\"Float64\" format=\"ascii\">\n");
  for (i = 0; i < s->nx; i++)
  {
    coor = s->xrange[0] + i * s->dx;
    fprintf(fp, "%1.6e ", coor);
  }
  fprintf(fp, "\n      </DataArray>\n");

  fprintf(fp, "      <DataArray Name=\"cy\" type=\"Float64\" format=\"ascii\">\n");
  for (i = 0; i < s->ny; i++)
  {
    coor = s->yrange[0] + i * s->dy;
    fprintf(fp, "%1.6e ", coor);
  }
  fprintf(fp, "\n      </DataArray>\n");

  fprintf(fp, "      <DataArray Name=\"cz\" type=\"Float64\" format=\"ascii\">\n");
  fprintf(fp, "0.0 ");
  fprintf(fp, "\n      </DataArray>\n");
  fprintf(fp, "      </Coordinates>\n");

  fprintf(fp, "      <CellData>\n");
  fprintf(fp, "      <DataArray Name=\"kappa_eff\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (i = 0; i < s->ncells; i++)
  {
    fprintf(fp, "%1.6e ", s->kappa[i]);
  }
  fprintf(fp, "\n      </DataArray>\n");
  fprintf(fp, "      </CellData>\n");

  fprintf(fp, "      <PointData>\n");
  fprintf(fp, "      <DataArray Name=\"h\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (i = 0; i < s->nnodes; i++)
  {
    fprintf(fp, "%1.6e ", s->h[i]);
  }
  fprintf(fp, "\n      </DataArray>\n");
  fprintf(fp, "      </PointData>\n");

  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </RectilinearGrid>\n");
  fprintf(fp, "</VTKFile>\n");

  fclose(fp);

  return (0);
}

ErrorCode FESLSPMCheckpoint(FESLSPM s)
{
  ErrorCode ierr;
  char name[128];

  sprintf(name, "spm-step%d.cpf", s->nt);
  ierr = FESLSPMView_Binary(s, name);

  return (0);
}

ErrorCode FESLSPMFieldView(FESLSPM s, FESLSPMFieldViewType type)
{
  ErrorCode ierr;
  char name[128];

  switch (type)
  {
  case FVT_ASCII:
    ierr = FESLSPMFieldView_ascii(s);
    break;
  case FVT_VTR:
    sprintf(name, "spm-step%d.vtr", s->nt);
    FESLSPMFieldView_VTR(s, name);
    break;
  case FVT_VTRB:
    sprintf(name, "spm-step%d.vtr", s->nt);
    FESLSPMFieldView_VTRB(s, name);
    break;
  default:
    sprintf(name, "spm-step%d.cpf", s->nt);
    ierr = FESLSPMView_Binary(s, name);
    break;
  }
  return (0);
}

ErrorCode FESLSPMLoadFromCheckpoint(FESLSPM s, const char name[])
{
  ErrorCode ierr;
  FILE *fp;

  fp = fopen(name, "r");
  if (!fp)
    ERRQ2("Unable to open file", name);

  fread(s->xrange, sizeof(double), 2, fp);
  fread(s->yrange, sizeof(double), 2, fp);

  fread(&s->time_end, sizeof(double), 1, fp);
  fread(&s->dt, sizeof(double), 1, fp);
  fread(&s->dt_const, sizeof(double), 1, fp);
  fread(&s->time, sizeof(double), 1, fp);

  fread(&s->nx, sizeof(int), 1, fp);
  fread(&s->ny, sizeof(int), 1, fp);
  fread(&s->mx, sizeof(int), 1, fp);
  fread(&s->my, sizeof(int), 1, fp);
  fread(&s->nnodes, sizeof(int), 1, fp);
  fread(&s->ncells, sizeof(int), 1, fp);

  fread(&s->dx, sizeof(double), 1, fp);
  fread(&s->dy, sizeof(double), 1, fp);
  fread(&s->dt_theta, sizeof(double), 1, fp);

  fread(&s->nt, sizeof(int), 1, fp);
  fread(&s->nt_io, sizeof(int), 1, fp);
  fread(&s->nt_max, sizeof(int), 1, fp);

  fread(&s->k0, sizeof(double), 1, fp);
  fread(&s->d_exp, sizeof(double), 1, fp);
  fread(&s->Sc, sizeof(double), 1, fp);

  {
    int dummy;

    fread(&dummy, sizeof(int), 1, fp);
    s->mode = (FESLSPMDiffusionMode)dummy;
    fread(&dummy, sizeof(int), 1, fp);
    s->dtmode = (FESLSPMTimestepMode)dummy;
  }

  s->h_ic = malloc(sizeof(double) * s->nnodes);
  s->h = malloc(sizeof(double) * s->nnodes);
  s->h_c = malloc(sizeof(double) * s->nnodes);
  s->u_v = malloc(sizeof(double) * s->nnodes);
  s->u_h = malloc(sizeof(double) * 2 * s->nnodes);
  s->k0_n = malloc(sizeof(double) * s->nnodes);
  s->d_exp_n = malloc(sizeof(double) * s->nnodes);
  s->Sc_n = malloc(sizeof(double) * s->nnodes);

  s->kappa = malloc(sizeof(double) * s->ncells);

  fread(s->h_ic, sizeof(double), s->nnodes, fp);
  fread(s->h, sizeof(double), s->nnodes, fp);
  fread(s->h_c, sizeof(double), s->nnodes, fp);
  fread(s->u_v, sizeof(double), s->nnodes, fp);
  fread(s->u_h, sizeof(double), 2 * s->nnodes, fp);
  fread(s->kappa, sizeof(double), s->ncells, fp);
  fread(s->k0_n, sizeof(double), s->nnodes, fp);
  fread(s->d_exp_n, sizeof(double), s->nnodes, fp);
  fread(s->Sc_n, sizeof(double), s->nnodes, fp);
  fread(s->bctype, sizeof(char), 4, fp);
  fread(&s->dt_min, sizeof(double), 1, fp);
  fread(&s->dt_max, sizeof(double), 1, fp);
  fclose(fp);

  return (0);
}

ErrorCode FESLSPMNodeDataGetRange(FESLSPM s, double x[], int start, int stride, double range[])
{
  int n, k;
  double val;

  range[0] = 1.0e32;
  range[1] = -1.0e32;

  n = s->nx * s->ny / stride;
  for (k = 0; k < n; k++)
  {
    val = x[start + k * stride];
    if (val < range[0])
      range[0] = val;
    if (val > range[1])
      range[1] = val;
  }

  return (0);
}

ErrorCode FESLSPMCellDataGetRange(FESLSPM s, double x[], int start, int stride, double range[])
{
  int n, k;
  double val;

  range[0] = 1.0e32;
  range[1] = -1.0e32;

  n = s->mx * s->my / stride;
  for (k = 0; k < n; k++)
  {
    val = x[start + k * stride];
    if (val < range[0])
      range[0] = val;
    if (val > range[1])
      range[1] = val;
  }

  return (0);
}

ErrorCode FESLSPMNodeDataSet(FESLSPM s, double x[], double val)
{
  int n, k;

  n = s->nx * s->ny;
  for (k = 0; k < n; k++)
  {
    x[k] = val;
  }
  return (0);
}

ErrorCode FESLSPMCellDataSet(FESLSPM s, double x[], double val)
{
  int n, k;

  n = s->mx * s->my;
  for (k = 0; k < n; k++)
  {
    x[k] = val;
  }
  return (0);
}

ErrorCode _FieldDataScale(FESLSPM s, int n, double x[], double val)
{
  int k;

  for (k = 0; k < n; k++)
  {
    x[k] = x[k] * val;
  }
  return (0);
}

ErrorCode _FieldDataShift(FESLSPM s, int n, double x[], double val)
{
  int k;

  for (k = 0; k < n; k++)
  {
    x[k] = x[k] + val;
  }
  return (0);
}

ErrorCode FESLSPMNodeDataScale(FESLSPM s, double x[], double val)
{
  int ierr;
  ierr = _FieldDataScale(s, s->nnodes, x, val);
  return (0);
}

ErrorCode FESLSPMCellDataScale(FESLSPM s, double x[], double val)
{
  int ierr;
  ierr = _FieldDataScale(s, s->ncells, x, val);
  return (0);
}

ErrorCode FESLSPMNodeDataShift(FESLSPM s, double x[], double val)
{
  int ierr;
  ierr = _FieldDataShift(s, s->nnodes, x, val);
  return (0);
}

ErrorCode FESLSPMCellDataShift(FESLSPM s, double x[], double val)
{
  int ierr;
  ierr = _FieldDataShift(s, s->ncells, x, val);
  return (0);
}

ErrorCode FESLSPMGetVerticalVelocity(FESLSPM s, int *m, int *n, double **data)
{
  if (m)
    *m = s->nx;
  if (n)
    *n = s->ny;
  if (data)
    *data = s->u_v;
  return (0);
}

ErrorCode FESLSPMGetHorizontalVelocity(FESLSPM s, int *m, int *n, double **data)
{
  if (m)
    *m = s->nx;
  if (n)
    *n = s->ny;
  if (data)
    *data = s->u_h;
  return (0);
}

ErrorCode FESLSPMGetKappa(FESLSPM s, int *m, int *n, double **data)
{
  if (m)
    *m = s->mx;
  if (n)
    *n = s->my;
  if (data)
    *data = s->kappa;
  return (0);
}

ErrorCode FESLSPMGetHeight(FESLSPM s, int *m, int *n, double **data)
{
  if (m)
    *m = s->nx;
  if (n)
    *n = s->ny;
  if (data)
    *data = s->h;
  return (0);
}

ErrorCode FESLSPMCreateNodalField(FESLSPM s, int ncomponents, double *nodedata[])
{
  double *data;

  if (!nodedata)
    ERRQ("A valid pointer must be provided for nodedata")

  data = malloc(sizeof(double) * s->nnodes * ncomponents);
  memset(data, 0, sizeof(double) * s->nnodes * ncomponents);

  *nodedata = data;
  return (0);
}

ErrorCode FESLSPMCreateCellField(FESLSPM s, int ncomponents, double *celldata[])
{
  double *data;

  if (!celldata)
    ERRQ("A valid pointer must be provided for celldata")

  data = malloc(sizeof(double) * s->ncells * ncomponents);
  memset(data, 0, sizeof(double) * s->ncells * ncomponents);

  *celldata = data;
  return (0);
}

#define SPMDestroyField(data) \
  {                           \
    free(*(data));            \
    *(data) = NULL;           \
  }

ErrorCode FESLSPMDestroyNodalField(FESLSPM s, double **data)
{
  if (!data)
    return (0);
  SPMDestroyField(data);
  return (0);
}

ErrorCode FESLSPMDestroyCellField(FESLSPM s, double **data)
{
  if (!data)
    return (0);
  SPMDestroyField(data);
  return (0);
}

ErrorCode FESLSPMCellGetNodeValues(FESLSPM s, int ci, int cj, int ncomponents, double nodalfield[], double vertexvalues[])
{
  int c;
  int cnodes[4];

  CellGetNodeIdx(s, ci, cj, cnodes);
  if (ncomponents == 1)
  {

    vertexvalues[0] = nodalfield[cnodes[0]];
    vertexvalues[1] = nodalfield[cnodes[1]];
    vertexvalues[2] = nodalfield[cnodes[2]];
    vertexvalues[3] = nodalfield[cnodes[3]];

    return (0);
  }

  for (c = 0; c < ncomponents; c++)
  {
    vertexvalues[ncomponents * 0 + c] = nodalfield[ncomponents * cnodes[0] + c];
    vertexvalues[ncomponents * 1 + c] = nodalfield[ncomponents * cnodes[1] + c];
    vertexvalues[ncomponents * 2 + c] = nodalfield[ncomponents * cnodes[2] + c];
    vertexvalues[ncomponents * 3 + c] = nodalfield[ncomponents * cnodes[3] + c];
  }
  return (0);
}

ErrorCode FESLSPMProjectNodalField(FESLSPM s, double nodedata[])
{
  double *denom, *numer;
  double cellnodevalues[4], navg;
  int n, k, cellnodes[4];
  int ierr;

  ierr = FESLSPMCreateNodalField(s, 1, &numer);
  ierr = FESLSPMCreateNodalField(s, 1, &denom);

  /* loop over cells */
  CentroidIteratorBegin(s)

      /* get nodal values for this cell */
      ierr = FESLSPMCellGetNodeValues(s, mi, mj, 1, nodedata, cellnodevalues);

  navg = 0.25 * (cellnodevalues[0] + cellnodevalues[1] + cellnodevalues[2] + cellnodevalues[3]);

  CellGetNodeIdx(s, mi, mj, cellnodes);

  /* accumulate N_i(xc)*phi_i and N_i(xc) */
  for (k = 0; k < 4; k++)
  {
    numer[cellnodes[k]] += navg;
    denom[cellnodes[k]] += 1.0;
  }

  CentroidIteratorEnd(s)

      for (n = 0; n < s->nnodes; n++)
  {
    nodedata[n] = numer[n] / denom[n];
  }

  ierr = FESLSPMDestroyNodalField(s, &denom);
  ierr = FESLSPMDestroyNodalField(s, &numer);

  return (0);
}

ErrorCode FESLSPMSetSemiLagrangeInterpolationType(FESLSPM spm, SLInterpolationType type)
{
  switch (type)
  {
  case SLINTERPOLATION_BILINEAR:
    printf("[FESLSPM] SemiLagrange interpolation: Linear\n");
    spm->semilagrange_interpolation_type = type;
    break;
  case SLINTERPOLATION_CSPLINE:
    printf("[FESLSPM] SemiLagrange interpolation: Cubic-spline <default>\n");
    spm->semilagrange_interpolation_type = type;
    break;
  case SLINTERPOLATION_CSPLINE_MONOTONE:
    printf("[FESLSPM] SemiLagrange interpolation: Cubic-spline with bound preserving limiter (quasi-monotone)\n");
    spm->semilagrange_interpolation_type = type;
    break;
  case SLINTERPOLATION_UPWIND_CSPLINE:
    printf("[FESLSPM] SemiLagrange interpolation: Cubic-spline with upwinding\n");
    spm->semilagrange_interpolation_type = type;
    break;
  case SLINTERPOLATION_UPWIND_CSPLINE_MONOTONE:
    printf("[FESLSPM] SemiLagrange interpolation: Cubic-spline with upwinding + bound preserving limiter\n");
    spm->semilagrange_interpolation_type = type;
    break;

  default:
    ERRQ("SemiLagrange interpolation: Not recognized. Valid types are (0 (linear); 1 (Cubic spline); 2 (Cubic spline + quasi-monotone limiter)");
    break;
  }

  return (0);
}
