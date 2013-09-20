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
      .def("GetConditionalP", &SSAgreement::GetConditionalP)
  ;

  register_ptr_to_python<SSAgreementPtr>();



}
