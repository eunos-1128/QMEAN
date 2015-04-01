#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <qmean/vec_list_magic.hh>

#include <qmean/distance_constraints.hh>


using namespace boost::python;





void export_DistanceConstraints()
{

  def("DistanceConstraints",&qmean::DistanceConstraints,(arg("alignment_handle")));
  def("SeqSim",&qmean::seq_sim,(arg("SImMat"),arg("Seq1"),arg("Seq2")));

  


  class_<qmean::IndexPair>("IndexPair",init<int,int>())
    .def_readwrite("index_one",&qmean::IndexPair::index_one)
    .def_readwrite("index_two",&qmean::IndexPair::index_two)
    .def(self_ns::str(self))
  ;





  class_<std::map<qmean::IndexPair,float > >("DistCMap")
    .def(map_indexing_suite<std::map<qmean::IndexPair,float> >())
  ;



}

