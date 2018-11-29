#include "CellPoissonHypreSolver.h"
#include "CellPoissonFACOps.h"
#include "SAMRAI/solv/FACPreconditioner.h"
#include "CellPoissonFACSolver.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>

using namespace SAMRAI;

class PfmFACSolver
{
public:
   PfmFACSolver(std::string name, const tbox::Dimension &dim,
                std::shared_ptr<tbox::Database> input_db);

   std::shared_ptr<CellPoissonFACSolver>& getCellPoissonFACSolver()
   {
      return d_fac_solver;
   }

private:

   std::shared_ptr<CellPoissonHypreSolver> d_hypre_poisson;
   std::shared_ptr<CellPoissonFACOps> d_fac_ops;
   std::shared_ptr<solv::FACPreconditioner> d_fac_precond;
   std::shared_ptr<CellPoissonFACSolver> d_fac_solver;
};

