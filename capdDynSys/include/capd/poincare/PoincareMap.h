/////////////////////////////////////////////////////////////////////////////
/// @file PoincareMap.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_POINCARE_MAP_H_
#define _CAPD_POINCARE_POINCARE_MAP_H_

#include <string>
#include "capd/dynset/AffineCoordinateChange.h"
#include "capd/dynset/C0Set.h"
#include "capd/dynset/C1Set.h"
#include "capd/dynset/C2Set.h"
#include "capd/dynset/CnSet.h"
#include "capd/poincare/PoincareException.h"
#include "capd/poincare/BasicPoincareMap.h"

namespace capd{

/// This namespace contains classes to compute Poincare Maps and Time Maps
namespace poincare{
/// @addtogroup poincare 
/// @{

/**
 *  PoicareMap class rigorously computes Poincare Map.
 *
 *  For given Dynamical System and section it rigorously computes first return map
 *  to section (Poincare Map) and its derivatives.
 *
 *  It uses BasicPoincareMap and adds error estimates to the results.
 *
 *  Let
 *    - \f$ \varphi(t, x): R \times R^n -> R^n \f$ be a dynamical system generated by given vector field.
 *    - \f$ s: R^n \to R \f$ be a section function,
 *    - \f$ S = \{x \in R^n : s(x) = 0\} \f$
 *    - \f$ P: S \to S \f$ be a Poincare Map
 *    - for given point \f$ x \in S \f$ let T(x) be first return time (in given crossing direction)
 *      i.e.  \f$ P(x) = \varphi(T(x), x) \in S \f$
 *
 *  In the following we denote by
 *    - dP the derivative of Poincare Map P  : \f$ dP(x) = \frac{\partial P(x)}{\partial x} \f$
 *    - dT the derivative of T(x) : \f$ dT(x)  = \frac{\partial T(x)}{\partial x} \f$
 *    - dF the derivative of the flow : \f$ dF(x) = \frac{\partial \varphi(T(x), x)}{\partial x} \f$
 *  Then
 *    \f$ dP = dF + \frac{\partial \varphi}{\partial t} dT \f$
 *
 *  Parameters:
 *   - DS    dynamical system
 *       - DS::MapType  vector field type
 *   - SectionT  scalar function describing section :
 */

template<typename SolverT, typename SectionT = AbstractSection<typename SolverT::MatrixType> >
class PoincareMap : public BasicPoincareMap<SolverT, SectionT>
{
public:
  typedef SolverT Solver;              ///< ODE solver
  typedef typename Solver::VectorFieldType VectorFieldType;             ///< vector field type
  typedef typename Solver::MatrixType MatrixType;
  typedef typename Solver::VectorType VectorType;
  typedef typename Solver::ScalarType ScalarType;
  typedef typename ScalarType::BoundType BoundType;      ///< a type of the end for intervals
  typedef typename Solver::SolutionCurve CurveType;
  typedef typename VectorType::size_type size_type;      ///< integral type used to index containers (vectors, matrices, etc)
  typedef dynset::AffineCoordinateChange<MatrixType> AffineCoordinateChange;
  typedef SectionT SectionType;    ///< type of function \f$ s: R^n \to R \f$ that represents Poincare section.
  typedef typename SectionType::SectionDerivativesEnclosureType SectionDerivativesEnclosureType;
  typedef capd::dynset::C0Set<MatrixType> C0Set;   ///< type of abstract base class for all C0 sets
  typedef capd::dynset::C1Set<MatrixType> C1Set;   ///< type of abstract base class for all C1 sets
  typedef capd::dynset::C2Set<MatrixType> C2Set;   ///< type of abstract base class for all C2 sets
  typedef capd::dynset::CnSet<MatrixType,0> CnSet; ///< type of abstract base class for all C2 sets
  typedef typename C2Set::HessianType HessianType; ///< type of data structure for storing of Hessians of multivariate mappings.
  typedef typename CnSet::JetType JetType;

  typedef typename BasicPoincareMap<Solver,SectionT>::CrossingDirection CrossingDirection; ///< section crossing direction; possible values are PlusMinus, Both(default) and MinusPlus.

/// constant timeStepCorrectionFactor is used as a correction factor to multiply
/// the time step in procedure reachSection when coming close to section
/// and actual time step is to large. The value 0.9 is just from numerical simulations
/// as a good working. It can be changed to a value from the interval (0,1)
  //static
  BoundType timeStepCorrectionFactor;
  //static
  int maxCorrectionAttempts;         ///< maximal number of correction attempts
  //static
  ScalarType minCrossingTimeStep;    ///< minimal time step during section crossing

  /// Constructs PoincareMap for given dynamical system and section
  PoincareMap( Solver & solver,
               SectionType & section,
               CrossingDirection direction = Both,
               const BoundType & errorTolerance = 0.1
    );

  /// Computes value of n-th iterate of the Poincare Map. The result is returned in given affine coordinate system.
  /// @param[in] theSet - set of initial conditions
  /// @param[in] x0 - origin of new coordinates
  /// @param[in] A - linear part of new coordinates
  /// @param[out] returnTime - bound for return time to the section
  /// @param[in] n - iterate of Poincare map to be computed
  /// @returns A*(P(theSet)-x0)
  template<typename T>
  VectorType operator()(T& theSet, const VectorType& x0, const MatrixType& A, ScalarType& out_returnTime, int n = 1);
  
  /// Computes value of n-th iterate of the Poincare Map. The result is returned in given affine coordinate system.
  /// @param[in] x - set of initial conditions
  /// @param[in] x0 - origin of new coordinates
  /// @param[in] A - linear part of new coordinates
  /// @param[out] returnTime - bound for return time to the section
  /// @param[in] n - iterate of Poincare map to be computed
  /// @returns A*(P(theSet)-x0)
  VectorType operator()(const VectorType& x, const VectorType& x0, const MatrixType& A, ScalarType& out_returnTime, int n = 1);
  
  
  /// Computes value of n-th iterate of the Poincare Map (used for C^0 sets) and returns bound of return time
  template<typename T>
  VectorType operator()(T& theSet, ScalarType& out_returnTime, int n = 1);

  /// Computes value of n-th iterate of the Poincare Map (used for C^0 sets) and returns bound of return time
  VectorType operator()(const VectorType& x, ScalarType& out_returnTime, int n = 1);

  
  /// Computes value of n-th iterate of the Poincare Map (used for C^0 sets)
  template<typename T>
  VectorType operator()(T& theSet, int n = 1);

  /// Computes value of n-th iterate of the Poincare Map
  VectorType operator()(VectorType x, int n = 1);

  
  /// Computes value of n-th iterate of Poincare Map and derivatives of the flow  (used for C^1 sets)
  template<typename T>
  VectorType operator()(T& theSet, MatrixType & dF, ScalarType& out_returnTime, int n = 1);

  /// Computes value of n-th iterate of Poincare Map and derivatives of the flow  (used for C^1 sets)
  VectorType operator()(const VectorType& x, MatrixType & dF, ScalarType& out_returnTime, int n = 1);

  
  /// Computes value of n-th iterate of Poincare Map and derivatives of the flow  (used for C^1 sets)
  template<typename T>
  VectorType operator()(T& theSet, MatrixType & dF, int n = 1);

  /// Computes value of n-th iterate of Poincare Map and derivatives of the flow  (used for C^1 sets)
  VectorType operator()(const VectorType& x, MatrixType & dF, int n = 1);

  
  /// Computes value of n-th iterate of Poincare Map, derivatives and hessian of the flow (used for C^2 sets)
  template<typename T>
  VectorType operator()(T& theSet, MatrixType& dF, HessianType& hessian, int n = 1);

  /// Computes value of n-th iterate of Poincare Map, derivatives and hessian of the flow (used for C^2 sets)
  template<typename T>
  VectorType operator()(T& theSet, MatrixType& dF, HessianType& hessian, ScalarType& out_returnTime, int n = 1);

  /// Computes n-th iterate of Poincare Map and derivatives of the flow to given order (used for C^n sets)
  template<typename T>
  VectorType operator()(T& theSet, typename T::JetType& result, int n = 1);

  /// Computes n-th iterate of Poincare Map and derivatives of the flow to given order (used for C^n sets)
  template<typename T>
  VectorType operator()(T& theSet, typename T::JetType& result, ScalarType& out_returnTime, int n = 1);

  const SectionDerivativesEnclosureType& getSectionDerivativesEnclosure() const {
    return sectionDerivativesEnclosure;
  }
  
protected:
  template<class T>
  VectorType computePoincareMap(T& originalSet, int n);

  /// Function checks if we crossed section and then returned in one step.
  /// In this case it throws an exception of PoincareException<T> type.
  template<typename T>
  void checkTransversability ( const T & theSet,  const VectorType & bounds);

  template<typename T>
  ScalarType getSign( const T & theSet,  VectorType & position, bool updatePosition, const VectorType & bounds);

  template<typename T>
  ScalarType getSign( const T & theSet,  VectorType & position, bool updatePosition);

  template<typename T>
  ScalarType getSign(const T & theSet);

  /// The function propagates the set of initial conditions \f$ theSet \f$ by the flow until its \f$n\f$-th intersection with the Poincare section.
  /// @return If succeed, \f$ theSet \f$ is one-integration-step close to the section and \f$ nextSet \f$ is either on the section or just after the section.
  /// @return true, if \f$ nextSet \f$ is after section and false, otherwise.
  /// @note: An exception is thrown if the solver cannot integrate the set \f$ theSet \f$ until its \f$ n \f$-th intersection with Poicnare section.
  template<typename T>
  void integrateUntilSectionCrossing(T& theSet, T& nextSet, int n = 1);

  /// On input: theSet is 'one step' before the section
  /// On output theSet is before the section and closer than sizeFactor*diam(theSet) in perpendicular direction to the section.
  template<typename T>
  void getCloserToSection(T& theSet, T& nextSet);

  /// Crosses the section and returns the value of Poincare Map.
  /// It updates also derivatives.
  /// We expect theSet to be just before the section.
  /// After this function theSet is on the section or just after the section.
  /// This function is common for all types of C^0, C^1, C^2 and C^n sets
  template<typename T>
  VectorType crossSection(T& theSet, T& nextSet);

  /// If it is possible to cross the section in one step, the function does this and returns bounf Computes return time by means of the Newton method. If not possible, returns false.
  template<typename T>
  bool crossSectionInOneStep(T& theSet, T& nextSet, ScalarType& oneStepReturnTime, VectorType& bound);

  MatrixType sectionCoordinateSystem;
  SectionDerivativesEnclosureType sectionDerivativesEnclosure;

  ScalarType m_signBeforeSection;
  ScalarType m_lastTimeStep;
}; // end of template PoincareMap

/// @}
}} // namespace capd::poincare

#include "capd/poincare/PoincareMap_templateMembers.h"
#include "capd/poincare/PoincareMap_templateOperator.h"

#endif  /* _POINCAREMAP_H */

