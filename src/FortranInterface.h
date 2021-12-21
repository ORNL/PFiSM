#ifndef PFiSM_FortranInterface_H
#define PFiSM_FortranInterface_H

#include "SAMRAI/SAMRAI_config.h"

extern "C" {
void SAMRAI_F77_FUNC(comprhs3d, COMPRHS3D)(const int&, const int&, const int&,
                                           const int&, const int&, const int&,
                                           const double*, const double*,
                                           const int&, const double&,
                                           const double*, const double&,
                                           const double&, double*);
void SAMRAI_F77_FUNC(comprhsphase3d, COMPRHSPHASE3D)(
    const int&, const int&, const int&, const int&, const int&, const int&,
    const double*, const int&, const double*, const int&, const double*,
    const double&, const double&, const double&, const double&, const double&,
    const double*);
void SAMRAI_F77_FUNC(compcforphase,
                     COMPCFORPHASE)(const int&, const int&, const int&,
                                    const int&, const int&, const int&,
                                    const double*, const int&, const double&,
                                    const double&, const double&, const double*,
                                    const int&);
void SAMRAI_F77_FUNC(addvel2rhs3d,
                     ADDVEL2RHS3D)(const int&, const int&, const int&,
                                   const int&, const int&, const int&,
                                   const double*, const double*, const int&,
                                   const double&, const double*);
void SAMRAI_F77_FUNC(comprhsex3d,
                     COMPRHSEX3D)(const int&, const int&, const int&,
                                  const int&, const int&, const int&,
                                  const double*, const double*, const int&,
                                  const double&, double*);
}
#endif
