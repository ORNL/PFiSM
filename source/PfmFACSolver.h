#include "SAMRAI/solv/CellPoissonHypreSolver.h"
#include "SAMRAI/solv/CellPoissonFACOps.h"
#include "SAMRAI/solv/FACPreconditioner.h"
#include "SAMRAI/solv/CellPoissonFACSolver.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"

#include <string>

using namespace SAMRAI;

class PfmFACSolver
{
public:
   PfmFACSolver(std::string name, const tbox::Dimension &dim,
                std::shared_ptr<tbox::Database> input_db);

   std::shared_ptr<solv::CellPoissonFACSolver>& getCellPoissonFACSolver()
   {
      return d_fac_solver;
   }

private:

   std::shared_ptr<solv::CellPoissonHypreSolver> d_hypre_poisson;
   std::shared_ptr<solv::CellPoissonFACOps> d_fac_ops;
   std::shared_ptr<solv::FACPreconditioner> d_fac_precond;
   std::shared_ptr<solv::CellPoissonFACSolver> d_fac_solver;
};

