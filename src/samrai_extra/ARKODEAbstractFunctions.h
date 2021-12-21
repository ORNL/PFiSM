/*************************************************************************
 * Inspired by ARKODEAbstractFunctions
 ************************************************************************/

#ifndef included_ARKODEAbstractFunctions
#define included_ARKODEAbstractFunctions

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/solv/SundialsAbstractVector.h"


/**
 * Class ARKODEAbstractFunctions is an abstract base class that defines
 * an interface for the user-supplied RHSFunction and preconditioner
 * routines to be used with ARKODE and CVSpgmr via the C++ wrapper
 * class ARKODESolver.  To use ARKODE with the C++ wrapper one must
 * derive a subclass of this base class and pass it into the ARKODESolver
 * constructor.  The pure virtual member functions in this interface are
 * used by ARKODE and CVSpgmr during the ODE integration process.  The
 * complete argument lists in the function signatures defined by ARKODE
 * for the user-supplied routines have been preserved for the most part.
 * In a few cases, some arguments do not appear in the function signatures
 * below since they are superfluous via this interface.
 *
 * @see ARKODESolver
 * @see SundialsAbstractVector
 */

class ARKODEAbstractFunctions
{
 public:
   /**
    * The constructor and destructor for ARKODEAbstractFunctions
    * is empty.
    */
   ARKODEAbstractFunctions();
   virtual ~ARKODEAbstractFunctions();

   /**
    * User-supplied right-hand side function evaluation.
    *
    * The function arguments are:
    *
    *
    *
    * - \b t        (INPUT) {current value of the independent variable}
    * - \b y        (INPUT) {current value of dependent variable vector}
    * - \b y_dot   (OUTPUT){current value of the derivative of y}
    *
    *
    *
    *
    * IMPORTANT: This function must not modify the vector y.
    */
   virtual int evaluateRHSFunction(
       double t, SAMRAI::solv::SundialsAbstractVector* y,
       SAMRAI::solv::SundialsAbstractVector* y_dot) = 0;

   /*
    * added Implicit
    */
   virtual int evaluateRHSFunctionImp(
       double t, SAMRAI::solv::SundialsAbstractVector* y,
       SAMRAI::solv::SundialsAbstractVector* y_dot) = 0;

   /*
    * added Explicit
    */
   virtual int evaluateRHSFunctionExp(
       double t, SAMRAI::solv::SundialsAbstractVector* y,
       SAMRAI::solv::SundialsAbstractVector* y_dot) = 0;

   /**
    * User-supplied function for setting up the preconditioner
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int CVSpgmrPrecondSet(double t,
                                 SAMRAI::solv::SundialsAbstractVector* y,
                                 SAMRAI::solv::SundialsAbstractVector* fy,
                                 int jok, int* jcurPtr, double gamma) = 0;

   /**
    * User-supplied function for setting up the preconditioner
    * to be used in the solution of the linear system that arises
    * during Newton iteration.
    */
   virtual int CVSpgmrPrecondSolve(double t,
                                   SAMRAI::solv::SundialsAbstractVector* y,
                                   SAMRAI::solv::SundialsAbstractVector* fy,
                                   SAMRAI::solv::SundialsAbstractVector* r,
                                   SAMRAI::solv::SundialsAbstractVector* z,
                                   double gamma, double delta, int lr) = 0;

   /*!
    * @brief User-supplied function to project the current solution and
    * estimated error onto the constraint manifold.
    *
    * The function arguments are:
    *
    * - \b t (INPUT) {current value of the independent variable}
    * - \b y (INPUT) {current value of the dependent variable vector}
    * - \b corr (OUTPUT) {correction such that y+corr is on the constraint
    *                     manifold}
    * - \b epsProj (INPUT) {WRMS norm tolerance for a nonlinear solver
    *                       iteration}
    * - \b err (INPUT/OUTPUT) {input: unprojected error , output: projected
    *                          error}
    *
    * IMPORTANT: This function must not modify the vector y.
    */
   virtual int applyProjection(double t,
                               SAMRAI::solv::SundialsAbstractVector* y,
                               SAMRAI::solv::SundialsAbstractVector* corr,
                               double epsProj,
                               SAMRAI::solv::SundialsAbstractVector* err) = 0;

   /*!
    * @brief User-supplied right-hand side function used in evaluating finite
    * difference Jacobian -vector products.
    *
    * The function arguments are:
    *
    * - \b t (INPUT) {current value of the independent variable}
    * - \b y (INPUT) {current value of the dependent variable vector}
    * - \b ydot (OUTPUT) {current value of the derivative of y}
    *
    * IMPORTANT: This function must not modify the vector y
    */
   virtual int evaluateJTimesRHSFunction(
       double t, SAMRAI::solv::SundialsAbstractVector* y,
       SAMRAI::solv::SundialsAbstractVector* y_dot) = 0;
};

#endif
