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

