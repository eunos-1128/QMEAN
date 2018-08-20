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

