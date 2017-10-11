#ifndef DisCo_HH
#define DisCo_HH

#include <vector>
#include <boost/shared_ptr.hpp>
#include <ost/seq/sequence_handle.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/geom/vec3.hh>
#include <ost/tri_matrix.hh>

namespace qmean{

class DisCoContainer;
typedef boost::shared_ptr<DisCoContainer> DisCoContainerPtr;


class DisCoContainer{
  
public:

  DisCoContainer(const ost::seq::SequenceHandle& seqres);

  ~DisCoContainer();

  void Save(const String& filename) const;

  static DisCoContainerPtr Load(const String& filename);

  void AddData(const ost::seq::AlignmentHandle& aln,
               const geom::Vec3List& positions,
               const std::vector<uint>& pos_seqres_mapping);

  void CalculateConstraints(Real dist_cutoff = 15.0, Real gamma = 70.0,
                            Real seqsim_clustering_cutoff = 0.4, 
                            Real bin_size = 1.0);

  bool HasConstraint(uint i, uint j) const;

  const std::vector<Real>& GetConstraint(uint i , uint j) const;

  std::vector<Real> GetScores(const ost::mol::EntityView& view) const;

  Real GetDistCutoff() const { return dist_cutoff_; }

  Real GetGamma() const { return gamma_; }

  Real GetSeqSimClusteringCutoff() const { return seqsim_clustering_cutoff_; }

  uint GetNumBins() const { return num_bins_; }

  Real GetBinSize() const { return bin_size_; }

private:

  bool readonly_;
  bool valid_constraints_;
  ost::seq::SequenceHandle seqres_;
  ost::seq::AlignmentList aln_;
  std::vector<geom::Vec3List> pos_;
  std::vector<std::vector<uint> > pos_seqres_mapping_;
  ost::TriMatrix<std::vector<Real>* > disco_scores_;
  // info on how the disco scores have been generated:
  Real dist_cutoff_;
  Real gamma_;
  Real seqsim_clustering_cutoff_;
  uint num_bins_;
  Real bin_size_;

};

} // ns

#endif
