// Copyright (c) 2013-2020, SIB - Swiss Institute of Bioinformatics and
// Biozentrum - University of Basel
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

