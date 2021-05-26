#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"

using namespace SAMRAI;

// Function adapted from SAMRAI::solv::CellPoissonFACOps::computeVectorWeights()

void computeVectorWeights(
   const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int weight_id,
   int coarsest_ln,
   int finest_ln)
{
   assert(hierarchy);
   assert(weight_id>=0);

   const tbox::Dimension dim(hierarchy->getDim());

   if (coarsest_ln == -1) coarsest_ln = 0;
   if (finest_ln == -1) finest_ln = hierarchy->getFinestLevelNumber();
   if (finest_ln < coarsest_ln) {
      TBOX_ERROR(
         "Illegal level number range.  finest_ln < coarsest_ln."
         << std::endl);
   }

   for (int ln = finest_ln; ln >= coarsest_ln; --ln) {

      /*
       * On every level, first assign cell volume to vector weight.
       */

      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator p(level->begin());
           p != level->end(); ++p) {
         const std::shared_ptr<hier::Patch>& patch = *p;
         std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
            SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
               patch->getPatchGeometry()));

         assert(patch_geometry);

         const double* dx = patch_geometry->getDx();
         double cell_vol = dx[0];
         if (dim > tbox::Dimension(1)) {
            cell_vol *= dx[1];
         }

         if (dim > tbox::Dimension(2)) {
            cell_vol *= dx[2];
         }

         std::shared_ptr<pdat::CellData<double> > w(
            SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
               patch->getPatchData(weight_id)));
         assert(w);
         w->fillAll(cell_vol);
      }

      /*
       * On all but the finest level, assign 0 to vector
       * weight to cells covered by finer cells.
       */

      if (ln < finest_ln) {

         /*
          * First get the boxes that describe index space of the next finer
          * level and coarsen them to describe corresponding index space
          * at this level.
          */

         std::shared_ptr<hier::PatchLevel> next_finer_level(
            hierarchy->getPatchLevel(ln + 1));
         hier::BoxContainer coarsened_boxes = next_finer_level->getBoxes();
         hier::IntVector coarsen_ratio(next_finer_level->getRatioToLevelZero());
         coarsen_ratio /= level->getRatioToLevelZero();
         coarsened_boxes.coarsen(coarsen_ratio);

         /*
          * Then set vector weight to 0 wherever there is
          * a nonempty intersection with the next finer level.
          * Note that all assignments are local.
          */

         for (hier::PatchLevel::iterator p(level->begin());
              p != level->end(); ++p) {

            const std::shared_ptr<hier::Patch>& patch = *p;
            for (hier::BoxContainer::iterator i = coarsened_boxes.begin();
                 i != coarsened_boxes.end(); ++i) {

               hier::Box intersection = *i * (patch->getBox());
               if (!intersection.empty()) {
                  std::shared_ptr<pdat::CellData<double> > w(
                     SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                        patch->getPatchData(weight_id)));
                  TBOX_ASSERT(w);
                  w->fillAll(0.0, intersection);

               }  // assignment only in non-empty intersection
            }  // loop over coarsened boxes from finer level
         }  // loop over patches in level
      }  // all levels except finest
   }  // loop over levels
}

