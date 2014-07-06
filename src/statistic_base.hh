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

using namespace ost;

namespace qmean{

class StatisticBase;
class StatisticContainer;
typedef boost::shared_ptr<StatisticBase> StatisticBasePtr;
typedef boost::shared_ptr<StatisticContainer> StatisticContainerPtr;

typedef enum  {
  Interaction, CBeta,
  Reduced, Packing, Torsion, CBPacking
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
  StatisticBase() {}

  virtual PotentialType GetType()=0;

  void Save(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds)=0;

};

}//namespace
#endif
