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

#ifndef POTENTIAL_BASE_HH
#define POTENTIAL_BASE_HH

#include <qmean/statistic_base.hh>
#include <qmean/module_config.hh>
#include <map>
#include <vector>
#include <qmean/float_histogram.hh>
#include <ost/io/binary_data_source.hh>
#include <ost/io/binary_data_sink.hh>
#include <ost/io/container_serialization.hh>
#include <boost/filesystem/convenience.hpp>
#include <boost/pointer_cast.hpp>
#include <ost/io/io_exception.hh>
#include <limits>
#include <stdint.h>
#include <boost/multi_array.hpp>



namespace qmean{

class PotentialBase;
class PotentialContainer;
typedef boost::shared_ptr<PotentialBase> PotentialBasePtr;
typedef boost::shared_ptr<PotentialContainer> PotentialContainerPtr;


class DLLEXPORT_QMEAN PotentialContainer{
public:
  PotentialContainer() { }

  static PotentialContainerPtr Load(const String& file_name);

  void Save(const String& file_name);

  void Remove(const String& key);

  std::vector<String> GetKeys();

  int GetSize() { return potentials_.size(); }

  PotentialBasePtr& operator [] (const String& key) { return potentials_[key]; }

private:
  std::map<String,PotentialBasePtr> potentials_;
};


class DLLEXPORT_QMEAN PotentialBase{
public:
  PotentialBase() { }

  virtual ~PotentialBase() = 0;

  virtual PotentialType GetType() const = 0;

  void Save(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds)=0;

};

}
#endif

