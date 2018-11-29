#include "PfmFACSolver.h"


PfmFACSolver::PfmFACSolver(std::string name,
                           const tbox::Dimension &dim,
                           std::shared_ptr<tbox::Database> input_db)
{
   std::shared_ptr<tbox::Database> obj_db(input_db->getDatabase(name));

   std::string hypre_poisson_name(name + "::hypre_solver");
   std::string fac_ops_name(name + "::fac_ops");
   std::string fac_precond_name(name + "::fac_precond");
   std::string fac_solver_name(name + "::fac_solver");

   d_hypre_poisson.reset( new CellPoissonHypreSolver(
            dim,
            hypre_poisson_name,
            obj_db->isDatabase("hypre_solver") ?
            obj_db->getDatabase("hypre_solver") :
            std::shared_ptr<tbox::Database>()) );

   d_fac_ops.reset( new CellPoissonFACOps(
            d_hypre_poisson,
            dim,
            fac_ops_name,
            obj_db->isDatabase("fac_ops") ?
            obj_db->getDatabase("fac_ops") :
            std::shared_ptr<tbox::Database>()) );

   d_fac_precond.reset( new solv::FACPreconditioner(
            fac_precond_name,
            d_fac_ops,
            obj_db->isDatabase("fac_precond") ?
            obj_db->getDatabase("fac_precond") :
            std::shared_ptr<tbox::Database>()) );
 
   d_fac_solver.reset( new CellPoissonFACSolver(
            dim,
            fac_solver_name,
            d_fac_precond,
            d_fac_ops,
            obj_db->isDatabase("fac_solver") ?
            obj_db->getDatabase("fac_solver") :
            std::shared_ptr<tbox::Database>()) );
}

