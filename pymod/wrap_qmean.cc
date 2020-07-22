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
#include <ost/mol/entity_handle.hh>
#include <qmean/version.hh>

void export_Torsion();
void export_Interaction();
void export_Packing();
void export_Reduced();
void export_CBeta();
void export_SSAgreement();
void export_Base();
void export_SphericalSmoother();
void export_CBPacking();
void export_HBond();
void export_disco();
void export_clash();
void export_gmqe();
void export_extract_data_helper();

using namespace boost::python;

BOOST_PYTHON_MODULE(_qmean)
{
  // attach version information to current modules scope
  scope().attr("qmean_version") = QMEAN_VERSION_STRING;
  scope().attr("qmean_version_major") = QMEAN_VERSION_MAJOR;
  scope().attr("qmean_version_minor") = QMEAN_VERSION_MINOR;
  scope().attr("qmean_version_patch") = QMEAN_VERSION_PATCH;

  // export qmean functionality
  export_Base();	
  export_Torsion();
  export_Interaction();
  export_Packing();
  export_Reduced();
  export_CBeta();
  export_SSAgreement();
  export_SphericalSmoother();
  export_CBPacking();
  export_HBond();
  export_disco();
  export_clash();
  export_gmqe();
  export_extract_data_helper();
}

