//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------


#include <boost/python.hpp>
#include <ost/db/linear_indexer.hh>
#include <ost/db/binary_container.hh>
#include <qmean/extract_data_helper.hh>

using namespace boost::python;
using namespace ost::db;


namespace{


template<typename T>
void VecToList(const std::vector<T>& v, list& l) {
  l = list();
  for(typename std::vector<T>::const_iterator i = v.begin(); 
  	  i != v.end(); ++i) {
  	l.append(*i);
  }
}

/////////////////////////////
// Wrapper for DisCo Stuff //
/////////////////////////////

// helper struct to reduce number of input parameters
struct DisCoDataContainers {

  DisCoDataContainers(const String& indexer_path,
                      const String& seqres_container_path,
                      const String& atomseq_container_path,
                      const String& position_container_path) {
    indexer = LinearIndexer::Load(indexer_path);
    seqres_container = LinearCharacterContainer::Load(seqres_container_path);
    atomseq_container = LinearCharacterContainer::Load(atomseq_container_path);
    position_container = LinearPositionContainer::Load(position_container_path);
  }

  LinearIndexerPtr indexer;
  LinearCharacterContainerPtr seqres_container;
  LinearCharacterContainerPtr atomseq_container;
  LinearPositionContainerPtr position_container;
};

// helper struct to wrap output

struct ExtractedDisCoData{
  std::vector<int> residue_numbers;
  geom::Vec3List ca_positions;
};

typedef boost::shared_ptr<ExtractedDisCoData> ExtractedDisCoDataPtr;
list WrapExtractedDisCoDataGetResNums(ExtractedDisCoDataPtr ptr) {
  list residue_numbers;
  VecToList(ptr->residue_numbers, residue_numbers); 
  return residue_numbers;
}

ExtractedDisCoDataPtr WrapExtractTemplateDataDisCo(const String& entry_name, 
                                                   const String& chain_name,
                                                   const ost::seq::AlignmentHandle& aln,
                                                   DisCoDataContainers& data_containers) {
  ExtractedDisCoDataPtr ptr(new ExtractedDisCoData);
  qmean::ExtractTemplateDataDisCo(entry_name, chain_name, aln, 
                                  data_containers.indexer, 
                                  data_containers.seqres_container, 
                                  data_containers.atomseq_container, 
                                  data_containers.position_container, 
                                  ptr->residue_numbers,
                                  ptr->ca_positions);
  return ptr;
} 

////////////////////////////
// Wrapper for GMQE Stuff //
////////////////////////////

// helper struct to reduce number of input parameters
struct GMQEDataContainers {

  GMQEDataContainers(const String& indexer_path,
                     const String& seqres_container_path,
                     const String& atomseq_container_path,
                     const String& dssp_container_path,
                     const String& n_position_container_path,
                     const String& ca_position_container_path,
                     const String& c_position_container_path,
                     const String& cb_position_container_path) {
    indexer = LinearIndexer::Load(indexer_path);
    seqres_container = LinearCharacterContainer::Load(seqres_container_path);
    atomseq_container = LinearCharacterContainer::Load(atomseq_container_path);
    dssp_container = LinearCharacterContainer::Load(dssp_container_path);
    n_position_container = LinearPositionContainer::Load(n_position_container_path);
    ca_position_container = LinearPositionContainer::Load(ca_position_container_path);
    c_position_container = LinearPositionContainer::Load(c_position_container_path);
    cb_position_container = LinearPositionContainer::Load(cb_position_container_path);
  }

  LinearIndexerPtr indexer;
  LinearCharacterContainerPtr seqres_container;
  LinearCharacterContainerPtr atomseq_container;
  LinearCharacterContainerPtr dssp_container;
  LinearPositionContainerPtr n_position_container;
  LinearPositionContainerPtr ca_position_container;
  LinearPositionContainerPtr c_position_container;
  LinearPositionContainerPtr cb_position_container;
};

// helper struct to wrap output
struct ExtractedGMQEData{
  std::vector<int> residue_numbers;
  String dssp;
  geom::Vec3List n_positions;
  geom::Vec3List ca_positions;
  geom::Vec3List c_positions;
  geom::Vec3List cb_positions; 
};
typedef boost::shared_ptr<ExtractedGMQEData> ExtractedGMQEDataPtr;
list WrapExtractedGMQEDataGetResNums(ExtractedGMQEDataPtr ptr) {
  list residue_numbers;
  VecToList(ptr->residue_numbers, residue_numbers);
  return residue_numbers;
}

ExtractedGMQEDataPtr WrapExtractTemplateDataGMQE(const String& entry_name, 
                                                 const String& chain_name,
                                                 const ost::seq::AlignmentHandle& aln,
                                                 GMQEDataContainers& data_containers) {
  ExtractedGMQEDataPtr ptr(new ExtractedGMQEData);
  qmean::ExtractTemplateDataGMQE(entry_name, chain_name, aln, 
                                 data_containers.indexer, 
                                 data_containers.seqres_container, 
                                 data_containers.atomseq_container,
                                 data_containers.dssp_container, 
                                 data_containers.n_position_container,
                                 data_containers.ca_position_container, 
                                 data_containers.c_position_container,
                                 data_containers.cb_position_container, 
                                 ptr->residue_numbers, ptr->dssp, 
                                 ptr->n_positions, ptr->ca_positions, 
                                 ptr->c_positions, ptr->cb_positions);
  return ptr;
}

}


void export_extract_data_helper() {

  class_<DisCoDataContainers>("DisCoDataContainers", init<const String&, 
                                                          const String&, 
                                                          const String&,
                                                          const String&>())
    .def_readonly("indexer", &DisCoDataContainers::indexer)
    .def_readonly("seqres_container", &DisCoDataContainers::seqres_container)
    .def_readonly("atomseq_container", &DisCoDataContainers::atomseq_container)
    .def_readonly("ca_position_container", &DisCoDataContainers::position_container)
  ;

  class_<ExtractedDisCoData, ExtractedDisCoDataPtr>("ExtractedDisCoData", no_init)
    .add_property("residue_numbers", &WrapExtractedDisCoDataGetResNums)
    .def_readonly("ca_positions", &ExtractedDisCoData::ca_positions)
  ;

  def("ExtractTemplateDataDisCo", &WrapExtractTemplateDataDisCo, (arg("entry_name"),
                                                                  arg("chain_name"),
                                                                  arg("aln"),
                                                                  arg("data_containers")));

  class_<GMQEDataContainers>("GMQEDataContainers", init<const String&, 
                                                        const String&, 
                                                        const String&,
                                                        const String&,
                                                        const String&,
                                                        const String&,
                                                        const String&,
                                                        const String&>())
    .def_readonly("indexer", &GMQEDataContainers::indexer)
    .def_readonly("seqres_container", &GMQEDataContainers::seqres_container)
    .def_readonly("atomseq_container", &GMQEDataContainers::atomseq_container)
    .def_readonly("dssp_container", &GMQEDataContainers::dssp_container)
    .def_readonly("n_position_container", &GMQEDataContainers::n_position_container)
    .def_readonly("ca_position_container", &GMQEDataContainers::ca_position_container)
    .def_readonly("c_position_container", &GMQEDataContainers::c_position_container)
    .def_readonly("cb_position_container", &GMQEDataContainers::cb_position_container)
  ;

  class_<ExtractedGMQEData, ExtractedGMQEDataPtr>("ExtractedGMQEData", no_init)
    .add_property("residue_numbers", &WrapExtractedGMQEDataGetResNums)
    .def_readonly("dssp", &ExtractedGMQEData::dssp)
    .def_readonly("n_positions", &ExtractedGMQEData::n_positions)
    .def_readonly("ca_positions", &ExtractedGMQEData::ca_positions)
    .def_readonly("c_positions", &ExtractedGMQEData::c_positions)
    .def_readonly("cb_positions", &ExtractedGMQEData::cb_positions)
  ;

  def("ExtractTemplateDataGMQE", &WrapExtractTemplateDataGMQE, (arg("entry_name"),
                                                                arg("chain_name"),
                                                                arg("aln"),
                                                                arg("data_containers")));
}

