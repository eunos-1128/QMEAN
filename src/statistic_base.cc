// QMEAN:
// ----------
// 
// Copyright (C) 2018 University of Basel and the Authors.
// Authors: Gabriel Studer, Pascal Benkert, Marco Biasini
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include <qmean/statistic_base.hh>

//shouldn't be done this way....
#include <qmean/reduced_statistics.hh>
#include <qmean/interaction_statistics.hh>
#include <qmean/cbeta_statistics.hh>
#include <qmean/packing_statistics.hh>
#include <qmean/torsion_statistics.hh>
#include <qmean/cb_packing_statistics.hh>
#include <qmean/hbond_statistics.hh>


using namespace ost;

namespace qmean{

StatisticContainerPtr StatisticContainer::Load(const String& filename){
  
  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open statistic container. File '"
       << filename << "' does not exist";
    throw io::IOException(ss.str());
  }
  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  StatisticContainerPtr p(new StatisticContainer());

  int size;
  String key;
  int type;

  ds & size;

  for(int i=0; i<size; ++i){
    ds & key;
    ds & type;

    switch(type){
      case Interaction:{
        InteractionStatisticPtr p_load(new InteractionStatistic);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      case CBeta:{
        CBetaStatisticPtr p_load(new CBetaStatistic);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      case Reduced: {
        ReducedStatisticPtr p_load(new ReducedStatistic);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      case Packing:{
        PackingStatisticPtr p_load(new PackingStatistic);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      case Torsion:{
        TorsionStatisticPtr p_load(new TorsionStatistic);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      case CBPacking:{
        CBPackingStatisticPtr p_load(new CBPackingStatistic);
        ds & *p_load;
        (*p)[key]=p_load;
        break;
      }
      case HBond:{
        HBondStatisticPtr p_load(new HBondStatistic);
        ds & *p_load;
        (*p)[key]=p_load;
      }
      default: throw io::IOException("Error in loading statistic container!");
    }
  }
  return p;
}

void StatisticContainer::Save(const String& file_name){

  std::ofstream stream(file_name.c_str(), std::ios_base::binary);
  std::stringstream ss;
  if(! stream){
    ss<<"Could not write to file '";
    ss<<file_name<<"'";
    throw ost::io::IOException(ss.str());
  }  
  io::BinaryDataSink ds(stream);

  int s = statistics_.size();
  ds & s;
  int t;

  for(std::map<String,StatisticBasePtr>::iterator it=statistics_.begin();it!=statistics_.end();++it){
    ds & it->first;
    t=it->second->GetType();
    ds & t;
    it->second->OnSave(ds);
  }
}

void StatisticContainer::Remove(const String& key){

  if(statistics_.count(key)>0){
    statistics_.erase(key);
    return;
  }
  LOG_ERROR("Key does not exist! Did not remove anything!");

}

std::vector<String> StatisticContainer::GetKeys(){
  std::vector<String> ret_vec;
  for(std::map<String,StatisticBasePtr>::iterator it=statistics_.begin();it!=statistics_.end();++it){
    ret_vec.push_back(it->first);
  }
  return ret_vec;
}

StatisticBase::~StatisticBase() { }

void StatisticBase::Save(const String& filename){
  std::ofstream stream(filename.c_str(), std::ios_base::binary);
  std::stringstream ss;
  if(! stream){
    ss<<"Could not write to file '";
    ss<<filename<<"'";
    throw ost::io::IOException(ss.str());
  }   
  io::BinaryDataSink ds(stream);
  this->OnSave(ds);
}

}//namespace

