#include "copyright.h"
/*============================================================================*/
/*! \file blast.c
 *  \brief Problem generator for spherical blast wave problem.
 *
 * PURPOSE: Problem generator for spherical blast wave problem.
 *
 * REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
 *   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.     */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:  */

/* This is not the default blast.c. I made some changes. */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  Prim1DS W;
  Cons1DS U1d;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real pressure,drat,prat,pa,da,x1,x2,x3;
  double theta;

  rin = 0.1;
  pa  = 0.1;
  da  = 1.0;
  drat = 1.0;
  prat = 1.1;

/* setup uniform ambient medium with spherical over-pressured region */

  W.d = da;
  W.Vx = 0.0;
  W.Vy = 0.0;
  W.Vz = 0.0;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
#ifndef ISOTHERMAL
        W.P = pa;
	if (x1 < 0.0) W.P = prat*pa;
#endif
        U1d = Prim1D_to_Cons1D(&(W),&Bx);

	pGrid->U[k][j][i].d  = U1d.d;
	pGrid->U[k][j][i].M1 = U1d.Mx;
	pGrid->U[k][j][i].M2 = U1d.My;
	pGrid->U[k][j][i].M3 = U1d.Mz;
#ifndef ISOTHERMAL
	pGrid->U[k][j][i].E  = U1d.E;
#endif
      }
    }
  }
#ifdef RESISTIVITY 
  eta_Ohm = 0.0;
  Q_AD    = par_getd("problem","Q_AD");
  Q_Hall  = 0.0;
  d_ind   = 0.0;
#endif


}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef RESISTIVITY
void get_eta_user(GridS *pG, int i, int j, int k,
                             Real *eta_O, Real *eta_H, Real *eta_A)
{

  *eta_O = 0.0;
  *eta_H = 0.0;
  *eta_A = 0.0;

  return;
}
#endif


void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}
