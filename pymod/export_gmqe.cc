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

#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include <qmean/gmqe_scores.hh>
#include <qmean/vec_list_magic.hh>


using namespace qmean;
using namespace boost::python;


namespace{

GMQEScoreCalculatorPtr WrapInit(PotentialContainerPtr potentials,
                                DisCoContainerPtr disco_container,
                                const ost::seq::SequenceHandle& seqres,
                                const ost::seq::SequenceHandle& psipred_pred,
                                const boost::python::list& psipred_cfi) {
  std::vector<int> v_psipred_cfi = ListToVec<int>(psipred_cfi);
  GMQEScoreCalculatorPtr p(new GMQEScoreCalculator(potentials, disco_container, 
                                                   seqres, psipred_pred,
                                                   v_psipred_cfi));
  return p;
}

boost::python::dict WrapEval(const GMQEScoreCalculatorPtr p,
                             const geom::Vec3List& n_positions,
                             const geom::Vec3List& ca_positions,
                             const geom::Vec3List& c_positions,
                             const geom::Vec3List& cb_positions,
                             const ost::seq::SequenceHandle& dssp_states,
                             const boost::python::list& residue_numbers) {
  std::vector<int> v_residue_numbers = ListToVec<int>(residue_numbers);
  Real dist_const, reduced, cb_packing, torsion, ss_agreement, coverage;
  p->Eval(n_positions, ca_positions, c_positions, cb_positions,
          dssp_states, v_residue_numbers, dist_const, reduced, cb_packing, 
          torsion, ss_agreement, coverage);
  boost::python::dict return_dict;
  return_dict["dist_const"] = dist_const;
  return_dict["reduced"] = reduced;
  return_dict["cb_packing"] = cb_packing;
  return_dict["torsion"] = torsion;
  return_dict["ss_agreement"] = ss_agreement;
  return_dict["coverage"] = coverage;

  return return_dict;
}



}


void export_gmqe() {

  class_<GMQEScoreCalculator>("GMQEScoreCalculator", no_init)
    .def("__init__", make_constructor(&WrapInit))
    .def("Eval", &WrapEval, (arg("n_positions"), arg("ca_positions"),
                             arg("c_positions"), arg("cb_positions"),
                             arg("dssp_states"), arg("residue_numbers")))
  ;

}
