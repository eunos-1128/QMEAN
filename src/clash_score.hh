#ifndef QMEAN_CLASH_SCORE_HH
#define QMEAN_CLASH_SCORE_HH

/*
  Author: Marco Biasini
 */

#include <ost/mol/entity_view.hh>
#include <ost/mol/entity_handle.hh>
#include <qmean/module_config.hh>

namespace qmean {

std::vector<Real> GetClashScores(const ost::mol::EntityView& target, 
                                 const ost::mol::EntityView& env);


} // ns

#endif

