#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>

#include <qmean/clash_score.hh>
#include <qmean/vec_list_magic.hh>

using namespace qmean;
using namespace boost::python;



namespace {
list WrapGetClashScores(ost::mol::EntityView& target, ost::mol::EntityView& env){
  std::vector<Real> scores = GetClashScores(target, env);
  list ret_list=VecToList<Real>(scores);
  return ret_list;
}
}

void export_clash() {
  def("GetClashScores", &GetClashScores, (arg("target"), arg("env")));
}
