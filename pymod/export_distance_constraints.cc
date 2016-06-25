#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <qmean/vec_list_magic.hh>
#include <qmean/distance_constraints.hh>


using namespace boost::python;

namespace{

  qmean::DCData WrapFillDCData(ost::seq::AlignmentHandle& msaln,boost::python::list& cluster, boost::python::list& filenames, ost::seq::alg::SubstWeightMatrixPtr subst, unsigned short dist_cutoff, unsigned short exp_factor){
    std::vector<std::vector<int> > v_cluster;
    boost::python::list cl;
    for(int i = 0; i < boost::python::len(cluster); ++i){      
      cl = boost::python::extract<boost::python::list>(cluster[i]);
      std::vector<int> v_cl = ListToVec<int>(cl);
      v_cluster.push_back(v_cl);
    }  
    std::vector<String> v_filenames = ListToVec<String>(filenames);
    
    return qmean::FillDCData(msaln,v_cluster, v_filenames, subst, dist_cutoff, exp_factor);
  }
}


void export_DistanceConstraints()
{
  def("FillDCData",&WrapFillDCData,(boost::python::arg("msaln"), boost::python::arg("cluster"), boost::python::arg("filenmes"),
                                    boost::python::arg("subst"), boost::python::arg("dist_cutoff")=15, boost::python::arg("exp_factor")=70));
  def("DCScore",&qmean::DCScore, (boost::python::arg("aln"), boost::python::arg("data"), boost::python::arg("dist_cutoff")=15,boost::python::arg("seq_sep")=0));
  def("SaveDCData", &qmean::SaveDCData);
  def("LoadDCData",&qmean::LoadDCData);
  def("DetermineFeatureValues",&qmean::DetermineFeatureValues);


  class_<qmean::DCData>("DCData", no_init)
    .def_readonly("seq",&qmean::DCData::seq)
    .def_readonly("cluster_seq_sim",&qmean::DCData::cluster_seq_sim)
    .def_readonly("cluster_seq_ids",&qmean::DCData::cluster_seq_ids)
    .def_readonly("cluster_weights",&qmean::DCData::cluster_weights)
    .def_readonly("dist_cutoff", &qmean::DCData::dist_cutoff)
    .def_readonly("exp_factor", &qmean::DCData::exp_factor);

  register_ptr_to_python<qmean::DCDataPtr>();
}

