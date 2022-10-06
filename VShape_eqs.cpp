
// #define XCat(tau, var) (Xhist->at_C(tau,var))
// #define XRat(tau, var) (Xhist->at_R(tau,var))

// #define XCat(tau, var) (Xhist->at_IntPol_3rd_C(tau,var))
// #define XRat(tau, var) (Xhist->at_IntPol_3rd_R(tau,var))

// #define XCat(tau, var) (Xhist->at_C_3rd_wSmallTau(tau,var))
// #define XRat(tau, var) (Xhist->at_R_3rd_wSmallTau(tau,var))

// #define XCat(tau, var) (Xhist->at_C_fast3rd(tau,var))
// #define XRat(tau, var) (Xhist->at_R_fast3rd(tau,var))

// #define XCat(tau, var) (Xhist->at_C_cubicHermite(tau,var))
// #define XRat(tau, var) (Xhist->at_R_cubicHermite(tau,var))

#define XCat(tau, var) (Xhist->at_C_fCH(tau,var))
#define XRat(tau, var) (Xhist->at_R_fCH(tau,var))


/////////////////////////shifted time
#ifdef SHIFTEDTIME
auto MLL_derivs = [](vars *X, vars_vec_wdX *Xhist, vars *d, parameters *p){
  d->G = -p->gamma_G * X->G + p->J_G - (sm::ExpR(X->G) - 1) * ( norm(X->E) + norm(XCat(2.0*p->T12,E)) * sm::ExpR(2.0*XRat(2.0*p->T12,Q)+XRat(2.0*p->T12,G)) );
  
//   d->Q = -p->gamma_Q * X->Q + p->J_Q - p->r_s * (sm::ExpR(2.0*X->Q) - 1) * sm::ExpR(X->G) * norm(X->E);
  d->Q = -p->gamma_Q * (X->Q - p->Q_0) - p->r_s * (sm::ExpR(2.0*X->Q) - 1) * sm::ExpR(X->G) * norm(X->E);
  
  d->E = -(p->delta + sm::img*p->omega)*X->E + p->delta * sqrt(p->kappa) * sm::ExpC( 0.5*(1.0-sm::img*p->alpha_G)*(XRat(2.0*p->T01,G)+XRat(p->T,G)) + (1.0-sm::img*p->alpha_Q)*XRat(p->T,Q) ) * XCat(p->T,E);
};
#endif

////////////////////////////real time
#ifdef REALTIME
auto MLL_derivs = [](vars *X, vars_vec_wdX *Xhist, vars *d, parameters *p){
  d->G = -p->gamma_G * X->G + p->J_G - (sm::ExpR(X->G) - 1) * ( norm( XCat(p->T01,E)) + norm(XCat(p->T01+2.0*p->T12,E)) * sm::ExpR(2.0*XRat(p->T12,Q)+XRat(2.0*p->T12,G)) );
  
//   d->Q = -p->gamma_Q * X->Q + p->J_Q - p->r_s * (sm::ExpR(2.0*X->Q) - 1) * sm::ExpR(XRat(p->T12,G)) * norm(XCat(p->T01+p->T12,E));
  d->Q = -p->gamma_Q * (X->Q - p->Q_0)- p->r_s * (sm::ExpR(2.0*X->Q) - 1) * sm::ExpR(XRat(p->T12,G)) * norm(XCat(p->T01+p->T12,E));
  
  d->E = -(p->delta + sm::img*p->omega)*X->E + p->delta * sqrt(p->kappa) * sm::ExpC( 0.5*(1.0-sm::img*p->alpha_G)*(XRat(p->T01,G)+XRat(p->T01+2.0*p->T12,G)) + (1.0-sm::img*p->alpha_Q)*XRat(p->T01+p->T12,Q) ) * XCat(p->T,E);
};
#endif








