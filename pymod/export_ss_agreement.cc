// Copyright (c) 2013-2018, SIB - Swiss Institute of Bioinformatics and
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

#include <qmean/ss_agreement.hh>

using namespace boost::python;
using namespace qmean;


void export_SSAgreement()
{
  class_<SSAgreement>("SSAgreement",init<>())

      .def(init<>())
      .def("GetPsiPredP", &SSAgreement::GetPsiPredP)
      .def("LogOddScore", &SSAgreement::LogOddScore)
      .def("GetDSSPP", &SSAgreement::GetDSSPP)
      .def("GetOverallP", &SSAgreement::GetOverallP)
  ;

  register_ptr_to_python<SSAgreementPtr>();
}

