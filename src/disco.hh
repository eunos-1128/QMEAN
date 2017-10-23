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


// struct to contain the constraint for a certain residue pair i,j
// plus some additional information
struct ConstraintInfo{

  Real max_seqsim;               // max avg sequence similarity of any of
                                 // the involved clusters
  Real max_seqid;                // max avg sequence identity of any of the
                                 // involved clusters
  Real variance;                 // variance calculated on ALL observed dists
  uint16_t num_clusters;         // number of involved clusters
  std::vector<Real> constraint;  // the actual constraint
};

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
                            Real seqsim_clustering_cutoff = 0.5, 
                            Real bin_size = 0.5);

  bool HasConstraint(uint i, uint j) const;

  const std::vector<Real>& GetConstraint(uint i , uint j) const;

  void FillData(const ost::mol::EntityView& view, 
                std::vector<Real>& scores) const;

  void FillData(const ost::mol::EntityView& view, 
                std::vector<Real>& scores,
                std::vector<uint>& counts,
                std::vector<Real>& avg_num_clusters,
                std::vector<Real>& avg_max_seqsim,
                std::vector<Real>& avg_max_seqid,
                std::vector<Real>& avg_variance) const;


  Real GetDistCutoff() const { return dist_cutoff_; }

  Real GetGamma() const { return gamma_; }

  Real GetSeqSimClusteringCutoff() const { return seqsim_clustering_cutoff_; }

  uint GetNumBins() const { return num_bins_; }

  Real GetBinSize() const { return bin_size_; }

private:

  void ClearConstraints();

  void FillViewData(const ost::mol::ResidueViewList& res_list,
                    geom::Vec3List& ca_positions,
                    std::vector<uint>& res_list_indices,
                    std::vector<uint>& seqres_indices) const;

  bool readonly_;
  bool valid_constraints_;
  ost::seq::SequenceHandle seqres_;
  ost::seq::AlignmentList aln_;
  std::vector<geom::Vec3List> pos_;
  std::vector<std::vector<uint> > pos_seqres_mapping_;
  ost::TriMatrix<ConstraintInfo*> constraint_infos_;
  // info on how the disco scores have been generated:
  Real dist_cutoff_;
  Real gamma_;
  Real seqsim_clustering_cutoff_;
  uint num_bins_;
  Real bin_size_;

};

} // ns

#endif
