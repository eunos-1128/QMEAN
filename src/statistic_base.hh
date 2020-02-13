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

#ifndef STATISTIC_BASE_HH
#define STATISTIC_BASE_HH


#include <qmean/module_config.hh>
#include <ost/io/binary_data_source.hh>
#include <ost/io/binary_data_sink.hh>
#include <ost/io/container_serialization.hh>
#include <boost/shared_ptr.hpp>
#include <ost/io/io_exception.hh>
#include <boost/filesystem/convenience.hpp>
#include <ost/log.hh>
#include <qmean/float_histogram.hh>
#include <vector>
#include <map>
#include <limits>
#include <exception>

using namespace ost;

namespace qmean{

class StatisticBase;
class StatisticContainer;
typedef boost::shared_ptr<StatisticBase> StatisticBasePtr;
typedef boost::shared_ptr<StatisticContainer> StatisticContainerPtr;

typedef enum  {
  Interaction, CBeta,
  Reduced, Packing, Torsion, CBPacking, HBond
} PotentialType;

class DLLEXPORT_QMEAN StatisticContainer{

public:

  StatisticContainer() { }

  static StatisticContainerPtr Load(const String& file_name);

  void Save(const String& file_name);

  void Remove(const String& key);

  std::vector<String> GetKeys();

  int GetSize() { return statistics_.size(); }

  StatisticBasePtr& operator[] (const String& key) { return statistics_[key]; }

private: 

  std::map<String,StatisticBasePtr> statistics_;
};


class DLLEXPORT_QMEAN StatisticBase{
public:
  StatisticBase() { }

  virtual ~StatisticBase() = 0;

  virtual PotentialType GetType() const = 0;

  void Save(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds) = 0;

};

}//namespace
#endif

