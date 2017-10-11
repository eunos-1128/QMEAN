#include <qmean/disco.hh>

#include <ost/seq/alg/merge_pairwise_alignments.hh>
#include <ost/seq/alg/subst_weight_matrix.hh>
#include <ost/mol/atom_view.hh>
#include <limits>
#include <queue>
#include <fstream>
#include <stdexcept>
#include <sstream>

namespace{

inline int ToIdx(char c) {
  if(c == '-') return 0;
  c = toupper(c);
  if(c < 'A' || c > 'Z') return 0;
  return static_cast<int>(toupper(c)-'A') + 1;
}


/* Takes a sequence alignment and a substitution matrix to calculate pairwise
sequence similarities. The first sequence in aln is assumed to be the SEQRES.
The pairwise sequence similarity among all sequences except SEQRES is filled
into the pairwise_seqsim variable. The sequence similarity of each sequence
relative to the SEQRES is filled into seqsim_to_seqres. */
void PairwiseSequenceSimilarities(const ost::seq::AlignmentHandle& aln,
                                  ost::seq::alg::SubstWeightMatrixPtr subst,
                                  ost::TriMatrix<Real>& pairwise_seqsim,
                                  std::vector<Real>& seqsim_to_seqres) {

  int num_sequences = aln.GetCount() - 1; // subtract SEQRES
  int aln_length = aln.GetLength();

  pairwise_seqsim = ost::TriMatrix<Real>(num_sequences, Real(0.0));

  // We directly transform the characters of the alignment into indices in
  // the substitution matrix (not including SEQRES).
  typedef std::vector<std::vector<char> > idx_mat_type; 
  idx_mat_type idx_mat(aln_length, std::vector<char>(num_sequences));
  for(int i = 0; i < num_sequences; ++i) {
    const String& s = aln.GetSequence(i + 1).GetString();
    for(int j = 0; j < aln_length; ++j) {
      idx_mat[j][i] = static_cast<char>(ToIdx(s[j]));
    }
  }

  // For faster lookup, we generate our own substitution matrix.
  // please note that we don't need any checks anymore for valid characters
  // in the sequences since we already transformed them to indices.
  // (everything invalid is simply 0, which corresponds to a gap)
  Real min_weight = subst->GetMinWeight();
  Real one_over_weight_range = subst->GetMaxWeight() - subst->GetMinWeight();
  one_over_weight_range = Real(1.0) / one_over_weight_range;

  // 'Z' - 'A' + 1 would be the size required to express the weight of all 
  // characters, we add one more to also cover '-'.
  const int own_subst_size = 'Z' - 'A' + 1 + 1;
  typedef std::vector<std::vector<int> > subst_mat_type;
  subst_mat_type own_subst(own_subst_size, std::vector<int>(own_subst_size, 0));
  int idx_one = 1;
  for(char one = 'A'; one <= 'Z'; ++one, ++idx_one) {
    int idx_two = 1;
    std::vector<int>& current_row = own_subst[idx_one];
    for(char two = 'A'; two <= 'Z'; ++two, ++idx_two) {
      current_row[idx_two] = static_cast<int>(subst->GetWeight(one, two));
    }
  }

  // Let's calculate the pairwise similarities!
  std::vector<int> summed_weights(num_sequences, 0);
  std::vector<int> counts(num_sequences, 0);
  for(int i = 0; i < num_sequences; ++i) {
    int num_other_sequences = num_sequences - i - 1;
    summed_weights.assign(num_other_sequences, 0);
    counts.assign(num_other_sequences, 0);
    // Instead of calculating the sequence similarity to all other sequences
    // separately, we process column by column. The advantage of that is that
    // we can fully skip the columns that contain a gap in the current sequence.
    for(int j = 0; j < aln_length; ++j) {
      const int ref_idx = idx_mat[j][i];
      if(ref_idx) {
        const char* idx_col = &(idx_mat[j][i+1]);
        const std::vector<int>& weight_col = own_subst[ref_idx];
        for(int k = 0; k < num_other_sequences; ++k) {
          summed_weights[k] += weight_col[idx_col[k]];
          counts[k] += (idx_col[k] > char(0)) ? 1 : 0;
        }
      }
    }

    for(int j = 0; j < num_other_sequences; ++j) {
      if(counts[j]){
        Real seq_sim = static_cast<Real>(summed_weights[j]) / counts[j];
        pairwise_seqsim.Set(i, i+1+j, 
                            (seq_sim - min_weight) * one_over_weight_range);
      }
    }
  }

  // lets fill in the sequence similarities to the SEQRES in a separate step
  // (code duplication I know)
  const String& seqres = aln.GetSequence(0).GetString();
  summed_weights.assign(num_sequences, 0);
  counts.assign(num_sequences, 0);
  seqsim_to_seqres.resize(num_sequences);

  for(int i = 0; i < aln_length; ++i) {
    int idx = ToIdx(seqres[i]);
    if(idx) {
      const char* idx_col = &(idx_mat[i][0]);
      const std::vector<int>& weight_col = own_subst[idx];
      for(int j = 0; j < num_sequences; ++j) {
        summed_weights[j] += weight_col[idx_col[j]];
        counts[j] += (idx_col[j] > char(0)) ? 1 : 0;
      }
    }
  }

  for(int i = 0; i < num_sequences; ++i) {
    if(counts[i]){
      Real seq_sim = static_cast<Real>(summed_weights[i]) / counts[i];
      seqsim_to_seqres[i] = (seq_sim - min_weight) * one_over_weight_range;
    }
  }
}


enum DistanceMetric {
  AVG_DIST_METRIC,
  MAX_DIST_METRIC,
  MIN_DIST_METRIC
};


struct Cluster {

  Cluster() { }

  void Merge(const Cluster& other){
    members.insert(members.end(), other.members.begin(), other.members.end());
  }

  std::vector<uint> members;
};


// gives avg distance between two clusters
Real ClusterDistanceAVG(const Cluster& one, 
                        const Cluster& two, 
                        const ost::TriMatrix<Real>& dist_matrix) {
  typedef std::vector<uint>::const_iterator MyIter;
  Real sum  = 0.0;
  for (MyIter i = one.members.begin(); i != one.members.end(); ++i) {
    for (MyIter j = two.members.begin(); j != two.members.end(); ++j) {
      sum += dist_matrix.Get(*i, *j);
    }
  }
  return sum / (one.members.size() * two.members.size());
}


// gives max distance between two clusters
Real ClusterDistanceMAX(const Cluster& one, 
                        const Cluster& two, 
                        const ost::TriMatrix<Real>& dist_matrix) {
  typedef std::vector<uint>::const_iterator MyIter;
  Real max_dist = -std::numeric_limits<Real>::max();
  Real current_dist;
  for (MyIter i = one.members.begin(); i != one.members.end(); ++i) {
    for (MyIter j = two.members.begin(); j != two.members.end(); ++j) {
      current_dist = dist_matrix.Get(*i, *j);
      if (current_dist > max_dist) max_dist = current_dist;
    }
  }
  return max_dist;
}


// gives min distance between two clusters
Real ClusterDistanceMIN(const Cluster& one, 
                        const Cluster& two, 
                        const ost::TriMatrix<Real>& dist_matrix) {
  typedef std::vector<uint>::const_iterator MyIter;
  Real min_dist = std::numeric_limits<Real>::max();
  Real current_dist;
  for (MyIter i = one.members.begin(); i != one.members.end(); ++i) {
    for (MyIter j = two.members.begin(); j != two.members.end(); ++j) {
      current_dist = dist_matrix.Get(*i, *j);
      if (current_dist < min_dist) min_dist = current_dist;
    }
  }
  return min_dist;
}


// returns whether actual val is smaller than extreme val
bool NegativeCompare(Real actual_val, Real extreme_val) {
  return actual_val < extreme_val;
}


// returns whether actual val is larger than extreme val
bool PositiveCompare(Real actual_val, Real extreme_val) {
  return actual_val > extreme_val;
}


struct QueueEntry{

  QueueEntry() { }
  QueueEntry(int a, int b, Real d): first(a), second(b), dist(d) { }

  int first;
  int second;
  Real dist;
};


// Similar compare functions as above but specifically for the priority queue
// The operators are the other way around than above! this is because of the 
// sorting behaviour of the queue
bool NegativeCompareQueueEntries(const QueueEntry& a, const QueueEntry& b){
  return a.dist > b.dist; 
}


bool PositiveCompareQueueEntries(const QueueEntry& a, const QueueEntry& b){
  return a.dist < b.dist; 
}


void DoClustering(const ost::TriMatrix<Real>& dist_matrix, 
                  DistanceMetric metric,
                  Real thresh, char direction, 
                  std::vector<std::vector<uint> >& output){

  output.clear();

  uint num = dist_matrix.GetSize();

  // do the stupid cases
  if(num == 0){
    return;
  }

  if(num == 1){
    std::vector<uint> super_big_cluster(1,0);
    output.push_back(super_big_cluster);
  }

  //let's get a function pointer to compare according to direction
  typedef bool (*cmp_f_ptr)(Real actual_val, Real other_val);
  cmp_f_ptr cmp_ptr;

  if (direction == '+') {
    cmp_ptr = &PositiveCompare;
  } else if (direction == '-') {
    cmp_ptr = &NegativeCompare;
  } else {
    throw std::runtime_error("Invalid direction observed in clustering!");
  }

  //let's get a function ptr to the right metric function
  typedef Real (*metric_f_ptr)(const Cluster& one, 
                               const Cluster& two, 
                               const ost::TriMatrix<Real>& dist_matrix);
  metric_f_ptr metric_ptr;

  switch (metric) {
    case AVG_DIST_METRIC: {
      metric_ptr = &ClusterDistanceAVG;
      break;
    }
    case MAX_DIST_METRIC: {
      metric_ptr = &ClusterDistanceMAX;
      break;
    }
    case MIN_DIST_METRIC: {
      metric_ptr = &ClusterDistanceMIN;
      break;
    }
    default: {
      metric_ptr = &ClusterDistanceAVG;
    }
  }

  // in agglomerative clustering we alway merge two clusters, we can therefore
  // represent this process using a binary tree, that has 2*num-1 nodes.
  // This is the upper limit for us because of the thresh value.
  // so we might stop merging before building up the whole tree
  // The idea is to build up sub trees from the bottom up and only
  // keep their roots (clusters that can be merged)
  std::vector<Cluster*> clusters(2*num - 1, NULL);
  for (uint i = 0; i < num; ++i) {
    clusters[i] = new Cluster;
    clusters[i]->members.push_back(i);
  }

  // generate a priority queue containing the distances between clusters
  typedef bool (*queue_cmp_f_ptr)(const QueueEntry&, const QueueEntry&);
  queue_cmp_f_ptr queue_cmp_ptr;
  if(direction == '+'){
    queue_cmp_ptr = &PositiveCompareQueueEntries;
  }
  else{
    queue_cmp_ptr = &NegativeCompareQueueEntries;
  } 
  std::priority_queue<QueueEntry,std::vector<QueueEntry>,queue_cmp_f_ptr> 
  distance_queue(queue_cmp_ptr);

  for(uint i = 0; i < num; ++i){
    for(uint j = i + 1; j < num; ++j){
      distance_queue.push(QueueEntry(i, j, dist_matrix.Get(i, j)));
    }
  }

  uint merge_iteration = 0;

  // there is a maximum of num - 1 merge events
  // after that, everything would be in the same cluster
  for( ; merge_iteration < num - 1; ++merge_iteration) {

    // search in priority queue until we get an element containing two
    // root clusters. Note, that we also find distance entries associated to 
    // clusters that have already been merged (and therefore deleted). 
    // We just ignore them and continue the search
    // Important: All clusters are merged into one cluster after num-1 merge steps. 
    // We're therefore guaranteed to find two clusters to merge at this point. 
    while((clusters[distance_queue.top().first] == NULL || 
           clusters[distance_queue.top().second] == NULL)){
      distance_queue.pop();
    }

    int cluster_a = distance_queue.top().first;
    int cluster_b = distance_queue.top().second;
    Real actual_distance = distance_queue.top().dist;
    distance_queue.pop();

    if(!(*cmp_ptr)(actual_distance, thresh)){
      // don't merge clusters anymore, they're all above/below thresh
      break;
    }

    // move cluster a to new location and merge cluster b into it
    int new_cluster_idx = num + merge_iteration;
    clusters[new_cluster_idx] = new Cluster(*clusters[cluster_a]);
    clusters[new_cluster_idx]->Merge(*clusters[cluster_b]);

    // delete the old clusters (we only keep the roots!)
    delete clusters[cluster_a];
    delete clusters[cluster_b];
    clusters[cluster_a] = NULL;
    clusters[cluster_b] = NULL;

    // calculate the distance from new cluster to all other root clusters
    for(int i = 0; i < new_cluster_idx; ++i){
      if(clusters[i] != NULL){
        Real distance = (*metric_ptr)(*clusters[i], *clusters[new_cluster_idx], 
                                      dist_matrix);
        distance_queue.push(QueueEntry(i, new_cluster_idx, distance));
      }
    }
  }

  //fill the output and clean up
  for (uint i = 0; i < num + merge_iteration; ++i) {
    if(clusters[i] != NULL){
      output.push_back(clusters[i]->members);
      delete clusters[i];
    }
  }
}


void EstimateClusterWeights(const std::vector<std::vector<uint> >& clusters, 
                            const std::vector<Real>& seqsim_to_target, 
                            Real gamma, std::vector<Real>& cluster_weights) {

  cluster_weights.assign(clusters.size(), Real(0.0));
  for(uint cluster_idx = 0; cluster_idx < clusters.size(); ++cluster_idx) {
    const std::vector<uint>& cluster = clusters[cluster_idx];
    uint cluster_size = cluster.size();
    if(cluster_size > 0) {
      for(std::vector<uint>::const_iterator it = cluster.begin(); 
          it != cluster.end(); ++it) {
        cluster_weights[cluster_idx] += seqsim_to_target[*it];
      }
      cluster_weights[cluster_idx] = 
      std::exp(gamma * cluster_weights[cluster_idx] / cluster_size);
    }
  }
}
  
void FillDistances(const ost::seq::SequenceHandle& seqres,
                   const std::vector<geom::Vec3List>& pos_list,
                   const std::vector<std::vector<uint> >& pos_seqres_mapping_list,
                   Real dist_cutoff,
                   ost::TriMatrix<std::vector<std::pair<uint16_t, uint16_t> >* >& 
                   pairwise_distances) {

  uint seqres_length = seqres.GetLength();
  pairwise_distances = 
  ost::TriMatrix<std::vector<std::pair<uint16_t, uint16_t> >* >(seqres_length,
                                                                NULL);

  Real squared_dist_cutoff = dist_cutoff * dist_cutoff;
  for(uint templ_idx = 0; templ_idx < pos_list.size(); ++templ_idx) {
    const geom::Vec3List& pos = pos_list[templ_idx];
    const std::vector<uint>& pos_mapping = pos_seqres_mapping_list[templ_idx];
    uint num_template_pos = pos.size();
    for(uint i = 0; i < num_template_pos; ++i) {
      for(uint j = i + 1; j < num_template_pos; ++j) {
        Real squared_dist = geom::Length2(pos[i] - pos[j]);
        if(squared_dist < squared_dist_cutoff) {
          uint r_idx_one = pos_mapping[i] - 1;
          uint r_idx_two = pos_mapping[j] - 1;
          std::vector<std::pair<uint16_t, uint16_t> >* dist_vec = 
          pairwise_distances.Get(r_idx_one, r_idx_two);
          if(dist_vec == NULL) {
            dist_vec = new std::vector<std::pair<uint16_t, uint16_t> >;
          } 
          uint16_t int_dist = 
          static_cast<uint16_t>(1000 * std::sqrt(squared_dist));
          dist_vec->push_back(std::make_pair(static_cast<uint16_t>(templ_idx),
                                             int_dist));
        }
      }
    }
  }
}


void FillDisCo(const std::vector<std::pair<uint16_t, uint16_t> >& dist,
               const std::vector<int>& cluster_assignments,
               const std::vector<Real>& cluster_weights,
               uint num_bins, Real bin_size, 
               std::vector<Real>& disco) {

  disco.assign(num_bins, Real(0.0));

  // count the number of occurences for each cluster
  std::vector<int> cluster_count(cluster_weights.size());
  for(uint i = 0; i < dist.size(); ++i) {
    uint16_t templ_idx = dist[i].first;
    cluster_count[cluster_assignments[templ_idx]] += 1;
  }

  // the per cluster weights of all clusters that occur at least
  // once must sum to one
  std::vector<Real> scaled_cluster_weights(cluster_weights.size());
  Real summed_cluster_weight = 0.0;
  for(uint i = 0; i < cluster_weights.size(); ++i) {
    if(cluster_count[i] > 0) {
      summed_cluster_weight += cluster_weights[i];
    }
  }

  for(uint i = 0; i < cluster_weights.size(); ++i) {
    scaled_cluster_weights[i] = cluster_weights[i] / summed_cluster_weight;
  }

  // for every observed distance we add a gaussian to the disco function
  // and estimate a prefactor that is determined by cluster size and 
  // scaled_cluster_weights
  for(uint i = 0; i < dist.size(); ++i) {

    uint cluster_idx = cluster_assignments[dist[i].first];
    Real d = static_cast<Real>(dist[i].second) * Real(0.001);
  
    // weighted by number of occurences in cluster
    Real weight = Real(1.0) / cluster_count[cluster_idx];
    // additional weight for according cluster
    weight *= scaled_cluster_weights[cluster_idx];

    for(uint j = 0; j < num_bins; ++j) {
      Real squared_diff = d - (j * bin_size);
      squared_diff *= squared_diff;
      disco[j] += weight * std::exp(Real(-0.5) * squared_diff); 
    }
  }
}


void FillDisCo(const ost::seq::SequenceHandle& seqres,
               const ost::TriMatrix<std::vector<std::pair<uint16_t, uint16_t> >* >& distances,
               const std::vector<int>& cluster_assignments,
               const std::vector<Real>& cluster_weights,
               uint num_bins, Real bin_size,
               ost::TriMatrix<std::vector<Real>* >& disco_scores) {

  uint num_templates = seqres.GetLength();
  for(uint i = 0; i < num_templates; ++i) {
    for(uint j = i + 1; j < num_templates; ++j) {
      const std::vector<std::pair<uint16_t, uint16_t> >* dist_vec = 
      distances.Get(i, j);
      if(dist_vec != NULL) {
        std::vector<Real>* disco_vec = disco_scores.Get(i,j);
        if(disco_vec == NULL) {
          disco_vec = new std::vector<Real>(num_bins, Real(0.0));
        }
        FillDisCo(*dist_vec, cluster_assignments, cluster_weights, 
                  num_bins, bin_size, *disco_vec);
        delete dist_vec;
      }
    }
  }
}


}

namespace qmean{


DisCoContainer::DisCoContainer(const ost::seq::SequenceHandle& seqres):
            readonly_(false),
            valid_constraints_(false),
            seqres_(seqres),
            disco_scores_(0, NULL),
            dist_cutoff_(std::numeric_limits<Real>::quiet_NaN()),
            gamma_(std::numeric_limits<Real>::quiet_NaN()),
            seqsim_clustering_cutoff_(std::numeric_limits<Real>::quiet_NaN()) {

  // check whether the seqres contains a gap
  if(seqres_.GetLength() != seqres_.GetGaplessString().size()) {
    throw std::runtime_error("The provided SEQRES must not contain a gap!");
  }

  // to save some memory we use datatypes, that cannot deal with
  // seqres_ sizes that cannot be represented by 16 bit integers
  // (which is super big anyway)
  if(seqres_.GetLength() >= std::numeric_limits<uint16_t>::max()) {
    throw std::runtime_error("Due to implementation reasons, the size of the "
                             "SEQRES must not exceed 65535!");
  } 
}


DisCoContainer::~DisCoContainer() {
  for(int i = 0; i < disco_scores_.GetSize(); ++i) {
    for(int j = i + 1; j < disco_scores_.GetSize(); ++j) {
      std::vector<Real>* vec = disco_scores_.Get(i, j);
      if(vec != NULL) {
        delete vec;
      }
    }
  }
}


void DisCoContainer::Save(const String& filename) const {

  if(!valid_constraints_) {
    throw std::runtime_error("There is nothing to save if you never called "
                             "the CalculateConstraints function!");
  }

  std::ofstream out_stream(filename.c_str(), std::ios::binary);
  if(!out_stream){
    throw std::runtime_error("Could not open file " + filename);
  }  

  // dump seqres
  uint16_t seqres_size = seqres_.GetLength();
  String seqres_string = seqres_.GetString();
  out_stream.write(reinterpret_cast<char*>(&seqres_size), sizeof(uint16_t));
  out_stream.write(&seqres_string[0], seqres_size);
  
  // dump parameters that have been used to calculate the DisCo scores
  float fvalue = dist_cutoff_;
  out_stream.write(reinterpret_cast<const char*>(&fvalue), sizeof(float));
  fvalue = gamma_;
  out_stream.write(reinterpret_cast<const char*>(&fvalue), sizeof(float));
  fvalue = seqsim_clustering_cutoff_;
  out_stream.write(reinterpret_cast<const char*>(&fvalue), sizeof(float));
  uint32_t ivalue = num_bins_;
  out_stream.write(reinterpret_cast<const char*>(&ivalue), sizeof(uint32_t));
  fvalue = bin_size_;
  out_stream.write(reinterpret_cast<const char*>(&fvalue), sizeof(float));

  // dump disco scores
  uint16_t N = disco_scores_.GetSize();
  uint num_bins = std::floor(dist_cutoff_ + 0.5) + 3;
  for(uint16_t i = 0; i < N; ++i) {
    for(uint16_t j = i + 1; j < N; ++j) {
      std::vector<Real>* vec = disco_scores_.Get(i, j);
      if(vec != NULL) {
        out_stream.write(reinterpret_cast<char*>(&i), sizeof(uint16_t));
        out_stream.write(reinterpret_cast<char*>(&j), sizeof(uint16_t));
        for(uint i = 0; i < vec->size(); ++i) {
          char s = std::floor((*vec)[i] * 100 + 0.5);
          out_stream.write(&s,1);
        }
      }
    }
  }
}


DisCoContainerPtr DisCoContainer::Load(const String& filename) {

  std::ifstream in_stream(filename.c_str(), std::ios::binary);
  if(!in_stream){
    throw std::runtime_error("Could not open file " + filename);
  }  

  // load the seqres
  uint16_t seqres_size;
  in_stream.read(reinterpret_cast<char*>(&seqres_size), sizeof(uint16_t));
  String s(seqres_size, '-');
  in_stream.read(&s[0], seqres_size);

  ost::seq::SequenceHandle seqres = ost::seq::CreateSequence("seqres", s);
  DisCoContainerPtr loaded_container(new DisCoContainer(seqres));
  loaded_container->readonly_ = true;

  // its only possible to save a DisCoContainer if the following flag is
  // true anyway
  loaded_container->valid_constraints_ = true;

  // load  parameters that have been used to calculate the DisCo scores
  Real fvalue;
  uint32_t ivalue;
  in_stream.read(reinterpret_cast<char*>(&fvalue), sizeof(float));
  loaded_container->dist_cutoff_ = fvalue;
  in_stream.read(reinterpret_cast<char*>(&fvalue), sizeof(float));
  loaded_container->gamma_ = fvalue;
  in_stream.read(reinterpret_cast<char*>(&fvalue), sizeof(float));
  loaded_container->seqsim_clustering_cutoff_ = fvalue;
  in_stream.read(reinterpret_cast<char*>(&ivalue), sizeof(uint32_t));
  loaded_container->num_bins_ = ivalue;
  in_stream.read(reinterpret_cast<char*>(&fvalue), sizeof(float));
  loaded_container->bin_size_ = fvalue;

  // load the disco scores
  loaded_container->disco_scores_ = 
  ost::TriMatrix<std::vector<Real>*>(seqres_size, NULL);
  
  // we just read and read until nothing comes anymore...
  uint16_t i;
  uint16_t j;
  String loaded_string(loaded_container->num_bins_, 'X');
  while(!in_stream.eof()) {
    in_stream.read(reinterpret_cast<char*>(&i), sizeof(uint16_t));
    in_stream.read(reinterpret_cast<char*>(&j), sizeof(uint16_t));
    in_stream.read(reinterpret_cast<char*>(&loaded_string[0]), 
                                           loaded_container->num_bins_);
    std::vector<Real>* loaded_constraint = 
    new std::vector<Real>(loaded_container->num_bins_);
    for(uint k = 0; k < loaded_container->num_bins_; ++k) {
      (*loaded_constraint)[k] = loaded_string[k] * 0.01;
    }
    loaded_container->disco_scores_.Set(i,j,loaded_constraint);
  }

  return loaded_container;
}

void DisCoContainer::AddData(const ost::seq::AlignmentHandle& aln,
                             const geom::Vec3List& positions,
                             const std::vector<uint>& pos_seqres_mapping) {

  if(readonly_) {
    throw std::runtime_error("Once you saved / loaded the DisCoContainer, "
                             "you cannot add more data!");
  }

  if(positions.size() != pos_seqres_mapping.size()) {
    throw std::runtime_error("the provided positions and the pos_seqres_mapping "
                             "must be consistent in size!");
  }

  int max_idx = seqres_.GetLength();
  for(uint i = 0; i < pos_seqres_mapping.size(); ++i) {
    if(pos_seqres_mapping[i] < 1 || pos_seqres_mapping[i] > max_idx) {
      std::stringstream ss;
      ss << "Entry with idx " << i << " in pos_seqres_mapping (";
      ss << pos_seqres_mapping[i]<<") does not fit into internal SEQRES!";
      throw std::runtime_error(ss.str());
    }
  }

  // to save some memory we use 16 bit integers as indices, check whether
  // we get bigger....
  if(pos_.size() >= std::numeric_limits<uint16_t>::max() - 1) {
    throw std::runtime_error("Due to implementation reasons, you can only add "
                             "a maximum of 65535 position vectors to "
                             "DisCoContainer");
  }

  aln_.push_back(aln);
  pos_.push_back(positions);
  pos_seqres_mapping_.push_back(pos_seqres_mapping);
}


void DisCoContainer::CalculateConstraints(Real dist_cutoff, Real gamma,
                                          Real seqsim_clustering_cutoff,
                                          Real bin_size) {

  if(dist_cutoff > 30.0) {
    // there is a limitation in one of the processing step, when representing
    // distances as integers to save some memory
    throw std::runtime_error("Due to implementation reasons, there is a "
                             "maximum of 30 Angstrom for the dist_cutoff!");
  }

  if(readonly_) {
    throw std::runtime_error("Once you saved / loaded the DisCoContainer, "
                             "you cannot recalculate the constraints anymore!");
  }

  dist_cutoff_ = dist_cutoff;
  gamma_ = gamma;
  seqsim_clustering_cutoff_ = seqsim_clustering_cutoff;
  bin_size_ = bin_size;

  // DisCo sums up gaussians with a standard deviation of one
  // to represent everything, we add a padding of three to allow
  // the constraint function to flatten out
  num_bins_ = std::ceil((dist_cutoff + 3.0) / bin_size_);

  ost::seq::AlignmentHandle full_aln = 
  ost::seq::alg::MergePairwiseAlignments(aln_, seqres_);

  // calculate pairwise sequence similarities
  ost::TriMatrix<Real> pairwise_seqsim(0);
  std::vector<Real> seqsim_to_target;
  ost::seq::alg::SubstWeightMatrixPtr subst(new ost::seq::alg::SubstWeightMatrix);
  subst->AssignPreset(ost::seq::alg::SubstWeightMatrix::BLOSUM62);
  PairwiseSequenceSimilarities(full_aln, subst, pairwise_seqsim, seqsim_to_target);

  // do the clustering
  std::vector<std::vector<uint> > clusters;
  DoClustering(pairwise_seqsim, AVG_DIST_METRIC, seqsim_clustering_cutoff_, 
               '+', clusters);
  std::vector<int> cluster_assignments(aln_.size());
  for(uint i = 0; i < clusters.size(); ++i) {
    for(uint j = 0; j < clusters[i].size(); ++j) {
      cluster_assignments[clusters[i][j]] = i;
    }
  }

  // estimate the per cluster weights
  std::vector<Real> cluster_weights;
  EstimateClusterWeights(clusters, seqsim_to_target, gamma_, cluster_weights);
  
  // estimate all the distances
  Real squared_dist_cutoff = dist_cutoff_ * dist_cutoff_;
  // the pair describes: (idx in pos_ / aln_, dist multiplied by 1000)
  ost::TriMatrix<std::vector<std::pair<uint16_t, uint16_t> >* > 
  pairwise_distances(0);
  FillDistances(seqres_, pos_, pos_seqres_mapping_, dist_cutoff_,
                pairwise_distances);

  // fill the disco scores
  // All the pointers to pairwise distances are deleted in 
  // this process... 
  // we just have to take care of the newly generated pointers in  disco_scores_ 
  // in the destructor of DisCoContainer
  FillDisCo(seqres_, pairwise_distances, cluster_assignments, cluster_weights,
            num_bins_, bin_size_, disco_scores_);

  valid_constraints_ = true;
}


bool DisCoContainer::HasConstraint(uint i, uint j) const {

  if(!valid_constraints_) {
    return false;
  }

  uint seqres_length = seqres_.GetLength();
  if(i < 1 || i > seqres_length || j < 1 || j > seqres_length) {
    throw std::runtime_error("Provided invalid resnum!");
  }

  return disco_scores_.Get(i-1,j-1) != NULL;
}


const std::vector<Real>& DisCoContainer::GetConstraint(uint i, uint j) const {

  if(!this->HasConstraint(i,j)) {
    throw std::runtime_error("No constraint observed for input resnums!");
  }

  std::vector<Real>* constraint = disco_scores_.Get(i-1, j-1);
  return *constraint;
}


std::vector<Real> DisCoContainer::GetScores(const ost::mol::EntityView& view) const {

  if(!valid_constraints_) {
    throw std::runtime_error("You must call the CalculateConstraints function "
                             "before you can get any scores!");
  }

  // only one chain allowed... the user has to apply a proper selection
  // statement in advance
  if(view.GetChainCount() != 1) {
    throw std::runtime_error("Input view must contain exactly 1 chain!");
  }

  // we extract the required data and at the same time check for 
  // consistency with the internal SEQRES
  ost::mol::ResidueViewList res_list = view.GetResidueList();
  geom::Vec3List ca_positions;
  std::vector<uint> res_list_indices;
  std::vector<uint> seqres_indices;
  uint seqres_length;

  for(uint i = 0; i < res_list.size(); ++i) {

    ost::mol::AtomView at = res_list[i].FindAtom("CA");

    if(at.IsValid()) {

      uint num = res_list[i].GetNumber().GetNum();
      if(num < 1 || num > seqres_length) {
        throw std::runtime_error("Observed a residue, that doesnt match the "
                                 "internal SEQRES!");
      }

      uint idx = num - 1;
      if(res_list[i].GetOneLetterCode() != seqres_[i]) {
        throw std::runtime_error("Observed a residue, that doesnt match the "
                                 "internal SEQRES!");
      }

      ca_positions.push_back(at.GetPos());
      res_list_indices.push_back(i);
      seqres_indices.push_back(num - 1);
    }
  }

  uint num_pos = ca_positions.size();
  std::vector<Real> summed_scores(num_pos, 0.0);
  std::vector<uint> counts(ca_positions.size(), 0);

  Real squared_dist_cutoff = dist_cutoff_ * dist_cutoff_;

  for(uint i = 0; i < num_pos; ++i) {
    for(uint j = i + 1; j < num_pos; ++j) {
      Real d = geom::Length2(ca_positions[i] - ca_positions[j]);
      if(d < squared_dist_cutoff) {
        
        std::vector<Real>* constraint = disco_scores_.Get(seqres_indices[i],
                                                          seqres_indices[j]);

        if(constraint != NULL) {

          d = std::sqrt(d);

          // we perform a linear interpolation
          uint bin_lower = static_cast<uint>(d / bin_size_);
          uint bin_upper = bin_lower + 1; 
          Real w_one = (bin_upper * bin_size_ - d) / bin_size_;
          Real w_two = Real(1.0) - w_one;
          Real score = w_one * (*constraint)[bin_lower] + 
                       w_two * (*constraint)[bin_upper];

          summed_scores[i] += score;
          summed_scores[j] += score;
          counts[i] += 1;
          counts[j] += 1;
        }
      }
    }
  }

  std::vector<Real> return_vec(res_list.size(), 
                               std::numeric_limits<Real>::quiet_NaN());

  for(uint i = 0; i < num_pos; ++i) {
    if(counts[i] > 0) {
      return_vec[res_list_indices[i]] = summed_scores[i] / counts[i];
    }
  }

  return return_vec;
}

} // ns
