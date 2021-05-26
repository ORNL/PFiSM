#include "HYPRE_struct_mv.h"


HYPRE_Int PFiSM_HYPRE_StructVectorGetBoxValues(HYPRE_StructVector vector,
                                               HYPRE_Int *ilower,
                                               HYPRE_Int *iupper,
                                               HYPRE_Complex *values);

HYPRE_Int PFiSM_HYPRE_StructVectorSetBoxValues(HYPRE_StructVector vector,
                                               HYPRE_Int *ilower,
                                               HYPRE_Int *iupper,
                                               HYPRE_Complex *values);

HYPRE_Int PFiSM_HYPRE_StructMatrixSetBoxValues(HYPRE_StructMatrix matrix,
                                               HYPRE_Int *ilower,
                                               HYPRE_Int *iupper,
                                               HYPRE_Int num_stencil_indices,
                                               HYPRE_Int *stencil_indices,
                                               HYPRE_Complex *values);
