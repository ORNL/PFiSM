#include "SAMRAI/hier/PatchHierarchy.h"

void computeVectorWeights(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int weight_id,
   int coarsest_ln,
   int finest_ln);

