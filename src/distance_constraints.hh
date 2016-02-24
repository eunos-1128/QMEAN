#include <ost/string_ref.hh>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <fstream>
#include <vector>
#include <ost/geom/vec3.hh>
#include <ost/geom/geom.hh>
#include <ost/mol/mol.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/subst_weight_matrix.hh>
#include <ost/seq/alg/sequence_similarity.hh>
#include <ost/seq/alg/sequence_identity.hh>
#include <math.h>


namespace qmean{

  struct DCData;
  typedef boost::shared_ptr<DCData> DCDataPtr;

  struct DCData{
    String seq;
    std::vector<float> cluster_weights;
    std::vector<float> cluster_seq_sim;
    std::vector<float> cluster_seq_ids;
    std::vector<std::vector<std::vector<float> > > tpl_distances;
    std::vector<std::vector<std::vector<unsigned short> > > cluster_sizes;
    unsigned short dist_cutoff;
    unsigned short exp_factor;
    // std::vector<std::vector<float> > max_scores;

    const std::vector<float>&  get_distances(int i, int j) const{
      if(i < j) return tpl_distances[i][j];
      return tpl_distances[j][i];
    };

    const std::vector<unsigned short>& get_cluster(int i, int j) const{
      if(i < j) return cluster_sizes[i][j];
      return cluster_sizes[j][i];
    }

    // const float& get_max_score(int i, int j) const{
    //   if(i < j) return max_scores[i][j];
    //   return max_scores[j][i];
    // }  
  };


  std::vector<std::vector<geom::Vec3> > ReadCalphaPos(const std::vector<String>& filenames,std::vector<bool>& valid_structures);
  std::vector<float> CalculateClusterSeqSim(const ost::seq::AlignmentHandle& msaln, const std::vector<std::vector<int> >& cluster,
                                            ost::seq::alg::SubstWeightMatrixPtr subst, const std::vector<bool>& valid_structures);
  std::vector<std::vector<std::vector<float> > > ExtractDistances(const ost::seq::AlignmentHandle& msaln,const std::vector<std::vector<int> >& cluster,
                                                                  const std::vector<std::vector<geom::Vec3> > tpl_ca_pos, const std::vector<bool>& valid_structures,
                                                                  std::vector<std::vector<std::vector<unsigned short> > >& cluster_sizes, unsigned short dist_cutoff);
  std::vector<float> DCScore(const ost::seq::AlignmentHandle& aln, const DCData& data, unsigned short dist_cutoff/*=15*/);
  void SaveDCData(DCData& data, const String& filename);
  DCData LoadDCData(const String& filename);
  DCData FillDCData(const ost::seq::AlignmentHandle& msaln,const std::vector<std::vector<int> >& cluster,const std::vector<String>& filenames,
                    ost::seq::alg::SubstWeightMatrixPtr subst, unsigned short dist_cutoff, unsigned short exp_factor);
  std::vector<std::vector<float> > DetermineFeatureValues(const ost::seq::AlignmentHandle& aln, const DCData& data);
  float CalculateVariance(const std::vector<float>& distances);
  float CalculateScoreIJ(float model_distance_ij ,const DCData& data, int i, int j);
  std::vector<float> CalculateClusterSeqIds(const ost::seq::AlignmentHandle& msaln, const std::vector<std::vector<int> >& cluster,
                                            const std::vector<bool>& valid_structures);
  std::vector<std::vector<int> > GetIndexMatrix(const ost::seq::AlignmentHandle& msaln);
  // std::vector<std::vector<float> > MaxScores(const DCData& data); 
}
