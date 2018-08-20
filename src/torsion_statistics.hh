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

#ifndef TORSION_STATISTICS_HH
#define TORSION_STATISTICS_HH

#include <qmean/statistic_base.hh>
#include <qmean/torsion_impl.hh>

namespace qmean {

class TorsionStatistic;
typedef boost::shared_ptr<TorsionStatistic> TorsionStatisticPtr;
typedef FloatHistogram<float,float,float,float,float,float> TorsionHistogram;

class DLLEXPORT_QMEAN TorsionStatistic : public StatisticBase, public impl::TorsionPotentialImpl {
public:
  TorsionStatistic() { }

  TorsionStatistic(std::vector<String>& group_identifier, std::vector<int>& num_of_bins);
 
  static TorsionStatisticPtr Load(const String& filename);

  virtual void OnSave(io::BinaryDataSink& ds);

  virtual void OnInteraction(const String& gi, std::vector<Real>& angles);

  void Extract(ost::mol::EntityView& target, Real weight);

  Real GetTotalCount() const;

  Real GetTotalCount(const String& gi) const;

  Real GetCount(const String& gi, std::vector<int>& bins) const;

  Real GetCount(std::vector<int>& bins) const;

  PotentialType GetType() const { return Torsion; }

  impl::TorsionOpts GetOpts() const { return opts_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & opts_;

    if(ds.IsSource()){
      for(std::vector<String>::iterator it=opts_.group_identifier.begin();it!=opts_.group_identifier.end();++it){
        histos_[*it]=TorsionHistogram();
      }
    }

    for(std::map<String,TorsionHistogram>::iterator it=histos_.begin();it!=histos_.end();++it){
      ds & it->second;
    }
  }

private:

  std::map<String, TorsionHistogram> histos_;
  Real weight_;
};

} //namespace


#endif // TORSION_STATISTICS_HH

