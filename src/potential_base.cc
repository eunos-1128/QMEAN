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

#include <qmean/potential_base.hh>

//I know, this is ugly as hell...
#include <qmean/reduced_potential.hh>
#include <qmean/interaction_potential.hh>
#include <qmean/cbeta_potential.hh>
#include <qmean/torsion_potential.hh>
#include <qmean/packing_potential.hh>
#include <qmean/cb_packing_potential.hh>
#include <qmean/hbond_potential.hh>

namespace qmean{

PotentialContainerPtr PotentialContainer::Load(const String& filename){
  
  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open potential container. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }
  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  PotentialContainerPtr p(new PotentialContainer);

  int size;
  String key;
  int type;

  ds & size;

  for(int i=0; i<size; ++i){
    ds & key;
    ds & type;

    switch(type){
      case Interaction:{
        InteractionPotentialPtr p_load(new InteractionPotential);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      case CBeta:{
        CBetaPotentialPtr p_load(new CBetaPotential);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      case Reduced: {
        ReducedPotentialPtr p_load(new ReducedPotential);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      case Packing:{
        PackingPotentialPtr p_load(new PackingPotential);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      case Torsion:{
        TorsionPotentialPtr p_load(new TorsionPotential);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      case CBPacking:{
        CBPackingPotentialPtr p_load(new CBPackingPotential);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      case HBond:{
        HBondPotentialPtr p_load(new HBondPotential);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      default: throw io::IOException("Error in loading potential container!");
    }
  }
  return p;
  
}


void PotentialContainer::Save(const String& filename){

  std::ofstream stream(filename.c_str(), std::ios_base::binary);
  std::stringstream ss;
  if(! stream){
    ss<<"Could not write to file '";
    ss<<filename<<"'";
    throw ost::io::IOException(ss.str());
  }  
  io::BinaryDataSink ds(stream);

  int s = potentials_.size();
  ds & s;
  int t;

  for(std::map<String,PotentialBasePtr>::iterator it=potentials_.begin();it!=potentials_.end();++it){
    ds & it->first;
    t=it->second->GetType();
    ds & t;
    it->second->OnSave(ds);
  }
  
}

void PotentialContainer::Remove(const String& key){

  if(potentials_.count(key)>0){
    potentials_.erase(key);
    return;
  }
  LOG_ERROR("Key does not exist! Did not remove anything!")

}

std::vector<String> PotentialContainer::GetKeys(){
  std::vector<String> ret_vec;
  for(std::map<String,PotentialBasePtr>::iterator it=potentials_.begin();it!=potentials_.end();++it){
    ret_vec.push_back(it->first);
  }
  return ret_vec;
}

PotentialBase::~PotentialBase() { }

void PotentialBase::Save(const String& filename){
 
  std::ofstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSink ds(stream);
  this->OnSave(ds);
}

}

