#ifndef PFiSM_Model_H
#define PFiSM_Model_H

#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/solv/SundialsAbstractVector.h"

using namespace SAMRAI;

class Model
{
 public:
   virtual void registerVisItDataWriter(
       std::shared_ptr<appu::VisItDataWriter>) = 0;
   virtual void setupSolutionVector(
       std::shared_ptr<hier::PatchHierarchy> hierarchy) = 0;
   virtual void setInitialConditions() = 0;
   virtual solv::SundialsAbstractVector* getSolutionVector() = 0;
   virtual double computeSolidFraction(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy) = 0;
   virtual void printCounters(const double) = 0;
};

#endif
