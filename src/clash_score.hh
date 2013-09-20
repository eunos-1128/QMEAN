#ifndef QMEAN_CLASH_SCORE_HH
#define QMEAN_CLASH_SCORE_HH

/*
  Author: Marco Biasini
 */

#include <ost/mol/entity_view.hh>
#include <ost/mol/entity_handle.hh>
#include <qmean/module_config.hh>

namespace qmean {

/// \defgroup Clash Steric Clash Score Calculation
/// 
/// The used scoring function is a simple steric energy function with linear 
/// fallof that is capped for small distances to avoid explosion of energies.
/// This is identical to the clash function used in SCWRL 3

//@{ 
/// \brief calculate clash score between two entity views
/// 
/// For each atom of ent_a the interaction with atoms of ent_b calculated.
/// \return 0.0 if there are no clashes and a positive clash score otherwise.
/// \sa \ref the_hammer.py "The Hammer Example"
Real DLLEXPORT_QMEAN ClashScore(const mol::EntityView& ent_a, const mol::EntityView& ent_b);

/// \brief calculate clash score between full entity and view
Real DLLEXPORT_QMEAN ClashScore(const mol::EntityHandle& ent_a, 
                                 const mol::EntityView& ent_b);
//// \brief calculate clash score of one single atom
/// 
/// \return floating point between 0 and 10
/// \sa \ref the_hammer.py "The Hammer Example"
Real DLLEXPORT_QMEAN ClashScore(const mol::AtomHandle& atom, 
                                 const mol::EntityView& ent_b);

/// \brief calculate steric energy of two atoms
Real DLLEXPORT_QMEAN StericEnergy(const geom::Vec3& pos1, Real r1,
                                   const geom::Vec3& pos2, Real r2);
//@}

/// \example the_hammer.py
/// 
/// Dynamic recalculation of clash score for a moving object. The real-valued
/// clash score is then color-mapped onto the objects. 
}

#endif
