#include <qmean/potential_base.hh>

//I know, this is ugly as hell...
#include <qmean/reduced_potential.hh>
#include <qmean/interaction_potential.hh>
#include <qmean/cbeta_potential.hh>
#include <qmean/torsion_potential.hh>
#include <qmean/packing_potential.hh>

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
      default: throw io::IOException("Error in loading potential container!");
    }
  }
  return p;
  
}


void PotentialContainer::Save(const String& filename){

  std::ofstream stream(filename.c_str(), std::ios_base::binary);
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

