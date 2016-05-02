#include <algorithm>
#include <qmean/distance_constraints.hh>
#include <ost/seq/aligned_column.hh>
#include <ost/seq/aligned_column_iterator.hh>
#include <numeric>
#include <boost/iostreams/filter/gzip.hpp> 
#include <boost/iostreams/filtering_stream.hpp>
#include <ost/io/io_exception.hh>


namespace qmean{
    std::vector<std::vector<geom::Vec3> > ReadCalphaPos(const std::vector<String>& filenames,
                                                        std::vector<bool>& valid_structures){

    std::vector<std::vector<geom::Vec3> > tpl_ca_pos; 
    for(int i=0; i<filenames.size(); ++i){
      std::vector<geom::Vec3>  ca_pos;
      String filename = filenames[i];
      String line;
      String old_res_num,res_num, atom_type;
      String x_str,y_str,z_str;
      Real x_pos,y_pos,z_pos;
      bool calpha_found = false;      
      bool all_ca_found = true;
      old_res_num = "none"; 

      boost::iostreams::filtering_stream<boost::iostreams::input> instream;
      instream.push(boost::iostreams::gzip_decompressor()); 
      const char * fn = filename.c_str();
      std::ifstream file(fn, std::ios_base::in | std::ios_base::binary);

      if(! file){
        std::stringstream ss;
        ss<<"Could not read file. File '";
        ss<<filename<<"' does not exist!";
        continue;
      }

      instream.push(file);
      while(std::getline(instream, line)) {
        if(line.size()>2 && line.substr(0,3) == "END" && calpha_found == false){all_ca_found = false;}
        if(line.size()>3 && line.substr(0,4) == "ATOM" || line.size()>5 && line.substr(0,6) == "HETATM"){       
          res_num = line.substr(22,4);
          atom_type = line.substr(12,4);
          res_num.erase(remove_if(res_num.begin(), res_num.end(), isspace),res_num.end());
          atom_type.erase(remove_if(atom_type.begin(), atom_type.end(), isspace),atom_type.end());

          if(atom_type == "CA"){
            x_str = line.substr(30,8);
            y_str = line.substr(38,8);
            z_str = line.substr(46,8);
            x_str.erase(remove_if(x_str.begin(), x_str.end(), isspace),x_str.end());
            y_str.erase(remove_if(y_str.begin(), y_str.end(), isspace),y_str.end());
            z_str.erase(remove_if(z_str.begin(), z_str.end(), isspace),z_str.end());
            x_pos = boost::lexical_cast<Real>(x_str);
            y_pos = boost::lexical_cast<Real>(y_str);
            z_pos = boost::lexical_cast<Real>(z_str);
            ca_pos.push_back(geom::Vec3(x_pos,y_pos,z_pos));
            calpha_found = true;
          }
          
          if(res_num != old_res_num && old_res_num != "none"){
            if(calpha_found == false){
              all_ca_found = false;
              // std::cout<<"ca missing "<<i<<" "<< filename<<std::endl;
              break;            
            }
            calpha_found = false;
          }
          old_res_num = res_num;
        }
      }
      file.close();
      if(all_ca_found){
        tpl_ca_pos.push_back(ca_pos);
        valid_structures.push_back(true);
      }
      else{
        tpl_ca_pos.push_back(std::vector<geom::Vec3>());
        valid_structures.push_back(false);
      }
    }    
    return tpl_ca_pos;
  }


  std::vector<float> CalculateClusterSeqSim(const ost::seq::AlignmentHandle& msaln, const std::vector<std::vector<int> >& cluster,
                                            ost::seq::alg::SubstWeightMatrixPtr subst, const std::vector<bool>& valid_structures){ 

    std::vector<float> cluster_similarities; //average sequence similarity of templates in clusters to target 
    float cl_sim;
    int tpl_idx;

    for(std::vector<std::vector<int> >::const_iterator it=cluster.begin();it!=cluster.end();++it){
      cl_sim = 0;
      for (std::vector<int>::const_iterator it2=(*it).begin();it2!=(*it).end();++it2){
        tpl_idx = *it2 +1;
        if(valid_structures[tpl_idx-1] == false) continue;
        cl_sim += ost::seq::alg::SequenceSimilarity(msaln,subst,true,0,tpl_idx);
      }
      cl_sim /= (*it).size();
      cluster_similarities.push_back(cl_sim);
    }
    return cluster_similarities;
  }

  
    std::vector<float> CalculateClusterSeqIds(const ost::seq::AlignmentHandle& msaln, const std::vector<std::vector<int> >& cluster,
                                            const std::vector<bool>& valid_structures){
    std::vector<float> cluster_seq_ids; //average sequnece identity of templates in clusters to target
    float cl_seq_id;
    int tpl_idx;

    for(std::vector<std::vector<int> >::const_iterator it=cluster.begin();it!=cluster.end();++it){
      cl_seq_id = 0;
      for (std::vector<int>::const_iterator it2=(*it).begin();it2!=(*it).end();++it2){
        tpl_idx = *it2 +1;
        if(valid_structures[tpl_idx-1] == false) continue;
        cl_seq_id += ost::seq::alg::SequenceIdentity(msaln,ost::seq::alg::RefMode::LONGER_SEQUENCE,0,tpl_idx);
      }
      cl_seq_id /= (*it).size();
      cluster_seq_ids.push_back(cl_seq_id);
    }
    return cluster_seq_ids;
  }

  
  std::vector<std::vector<int> > GetIndexMatrix(const ost::seq::AlignmentHandle& msaln){

    std::vector<std::vector<int> > idxs(msaln.GetCount());
    int c, i = 0;
    for( std::vector<std::vector<int> >::iterator it = idxs.begin(); it != idxs.end(); ++it, ++i){
      it->resize(msaln.GetLength());
      int offset = msaln.GetSequence(i).GetOffset();
      c = 0;
      ost::seq::ConstSequenceHandle seq = msaln.GetSequence(i);
      for(int j = 0; j<seq.GetLength(); j++){
        if(seq[j] != '-'){
          // for target sequence
          if(i == 0) (*it)[j] = c;
          // for templates the offset needs to be added
          else (*it)[j] = c+offset;
          c++;
        }
      }
    }
    return idxs;
  }


  std::vector<std::vector<std::vector<float> > > ExtractDistances(const ost::seq::AlignmentHandle& msaln,const std::vector<std::vector<int> >& cluster,   const std::vector<std::vector<geom::Vec3> > tpl_ca_pos,
                                                                  const std::vector<bool>& valid_structures, std::vector<std::vector<std::vector<unsigned short> > >& cluster_sizes, unsigned short dist_cutoff){

    std::vector<float> tpl_distances_ij;           // distances between reidue pair (i,j) in all templates
    std::vector<unsigned short> cluster_sizes_ij;  // info about number of distances in clusters for residue pair (i,j)
    float dist_ijk;                                // distance between reisdue pair (i,j) of template k
    int idx_i, idx_j, i = 0, j, k;
    geom::Vec3 ca_pos_i, ca_pos_j;

    // create a matrix to be able to look up c-alpha positios of each residue in the alignemnt
    std::vector<std::vector<int> > idxs = GetIndexMatrix(msaln);

    int n = msaln.GetSequence(0).GetGaplessString().size();
    std::vector<std::vector<std::vector<float> > > tpl_distances(n);
    for(std::vector<std::vector<std::vector<float> > >::iterator it = tpl_distances.begin(); it != tpl_distances.end(); ++it){
      it->resize(n);
    } 
    cluster_sizes.resize(n);
    for(std::vector<std::vector<std::vector<unsigned short> > >::iterator it = cluster_sizes.begin(); it != cluster_sizes.end(); ++it){
      it->resize(n);
    } 

    for (ost::seq::AlignmentHandle::iterator aln_it1 = msaln.begin(), e1 = msaln.end(); aln_it1 != e1; aln_it1++, i++) {
      const ost::seq::AlignedColumn& col_i = *aln_it1;
      j = i+1;

      if (col_i[0] == '-') continue;      
      for (ost::seq::AlignmentHandle::iterator aln_it2 = ost::seq::AlignedColumnIterator(msaln,j,msaln.GetLength()), e2 = msaln.end(); aln_it2 != e2; aln_it2++, j++) {            
        const ost::seq::AlignedColumn& col_j = *aln_it2;

        if (col_j[0] == '-') continue;
        tpl_distances_ij.clear();
        cluster_sizes_ij.clear();            
        for(std::vector<std::vector<int> >::const_iterator cl_it1=cluster.begin();cl_it1!=cluster.end();++cl_it1){
          for (std::vector<int>::const_iterator cl_it2=(*cl_it1).begin();cl_it2!=(*cl_it1).end();++cl_it2){
            k = (*cl_it2) ;                
            if ((col_i[k+1] == '-') || (col_j[k+1] == '-') || valid_structures[k] == false) continue;
            idx_i = idxs[k+1][i]; 
            idx_j = idxs[k+1][j];
            ca_pos_i = tpl_ca_pos[k][idx_i];
            ca_pos_j = tpl_ca_pos[k][idx_j];
            dist_ijk = geom::Distance(ca_pos_i,ca_pos_j);   

            if(dist_ijk < dist_cutoff){         
              tpl_distances_ij.push_back(dist_ijk);
            }                   
          }
          cluster_sizes_ij.push_back(tpl_distances_ij.size());
        }

        idx_i = idxs[0][i]; 
        idx_j = idxs[0][j];
        if(tpl_distances_ij.size()==0 ) continue;
        tpl_distances[idx_i][idx_j] = tpl_distances_ij;
        cluster_sizes[idx_i][idx_j] = cluster_sizes_ij;
      }
    }
    return tpl_distances;
  }


  float CalculateScoreIJ(float model_distance_ij ,const DCData& data, int i, int j){

    const float con = 0.4;  //for sigma = 1;
    std::vector<float> gaussians; 
    std::vector<float> tpl_distances_ij = data.get_distances(i,j);
    std::vector<unsigned short> cluster_sizes = data.get_cluster(i,j);      
    for(int k = 0; k< tpl_distances_ij.size(); k++){
      gaussians.push_back(con*exp(-pow(model_distance_ij-tpl_distances_ij[k],2)/2.));
    }
    float score_ij = 0;
    float cluster_weight_tot = 0;
    int cluster_start = 0;

    for(int k = 0; k< cluster_sizes.size(); k++){
      int cluster_end = cluster_sizes[k];
      if(cluster_start == cluster_end) continue;

      // sum up all gaussians in a cluster
      float cluster_gauss_sum = std::accumulate(gaussians.begin()+(cluster_start), gaussians.begin()+(cluster_end), 0.0);

      // divide by number of occurences of res pair in cluster
      cluster_gauss_sum /= float(cluster_end-cluster_start);  

      // weight with sequence similarity dependent weighting factor
      float cluster_weight = data.cluster_weights[k];
      cluster_gauss_sum *= cluster_weight;
      cluster_weight_tot += cluster_weight;

      // add cluster score to total score of residue pair (i,j)
      score_ij += cluster_gauss_sum;

      cluster_start = cluster_end;
    }
    // weights of clusters have to add up to one
    score_ij /= cluster_weight_tot;  
  return score_ij; 
  }    


  std::vector<float> DCScore(const ost::seq::AlignmentHandle& aln, const DCData& data ,unsigned short dist_cutoff/*=15*/){

    int j, i =0;
    float score_ij, model_distance_ij;
    std::vector<float> model_scores;

    float scores[aln.GetLength()]; 
    int num_inter[aln.GetLength()];
    memset( scores, 0, aln.GetLength()*sizeof(float));  
    memset( num_inter, 0, aln.GetLength()*sizeof(int)); 

    // check for valid DCData amd valid alignment
    if(aln.GetSequence(0).GetString()!= data.seq.c_str() ){
      throw std::invalid_argument("SEQRES of DCData does not match SEQRES in alignment!");
    }

    if(not aln.GetCount()>=2){
      throw std::invalid_argument("Number of sequences in alignment has to be at least 2!");
    }

    //Compute DisCo score
    for (ost::seq::AlignmentHandle::iterator aln_it1 = aln.begin(), e1 = aln.end(); aln_it1 != e1; aln_it1++, i++) {
      const ost::seq::AlignedColumn& col_i = *aln_it1;
      j = i+1;

      if (col_i[1] == '-') continue;
      for (ost::seq::AlignmentHandle::iterator aln_it2 = ost::seq::AlignedColumnIterator(aln,j,aln.GetLength()), e2 = aln.end(); aln_it2 != e2; aln_it2++, j++) {            
        const ost::seq::AlignedColumn& col_j = *aln_it2;
        if(col_j[1] == '-') continue;
        // if(std::abs(i-j)<=1) continue;

        ost::mol::ResidueView res_i = aln.GetSequence(1).GetResidue(i);
        ost::mol::AtomView ca_i = res_i.FindAtom("CA");
        ost::mol::ResidueView res_j = aln.GetSequence(1).GetResidue(j);             
        ost::mol::AtomView ca_j = res_j.FindAtom("CA"); 
        if(!ca_i.IsValid() || !ca_j.IsValid()) continue; 
        model_distance_ij = geom::Distance(ca_i.GetPos(),ca_j.GetPos());
        if(model_distance_ij >= dist_cutoff) continue;
        if(data.get_distances(i,j).empty()) continue;  //no templates for this residue pair

        float score_ij = CalculateScoreIJ(model_distance_ij, data, i, j);
        scores[i] += score_ij;
        scores[j] += score_ij;
        num_inter[i] += 1;
        num_inter[j] += 1;
      }
    }

    int first_res = aln.GetSequence(1).GetFirstNonGap();
    int last_res = aln.GetSequence(1).GetLastNonGap();
    
    for(int k = first_res; k<=last_res; k++ ){
      if(aln.GetSequence(1)[k] == '-') continue;
      if(num_inter[k] == 0)  model_scores.push_back(NAN);
      else{model_scores.push_back(scores[k]/(num_inter[k]*0.4));}
    }
  return model_scores;
  }


  // std::vector<std::vector<float> > MaxScores(const DCData& data){ 
  //
  //   int n = data.tpl_distances.size();
  //   std::vector<std::vector<float> > max_scores(n);
  //   for(std::vector<std::vector<float> >::iterator it = max_scores.begin(); it != max_scores.end(); ++it){
  //     it->resize(n);
  //   } 
  //   // find maximal possible score for residue pair ij  
  //   for(int i =0; i< n; i++){
  //     for(int j =i+1; j< n; j++){
  //       if(data.get_distances(i,j).empty()) continue;             
  //       float max_score_ij = 0;
  //       float stepsize = 0.1;
  //       std::vector<float> distances = data.get_distances(i,j);
  //       float min = *std::min_element(distances.begin(),distances.end());
  //       float max = *std::max_element(distances.begin(),distances.end());
  //       if(max == min){
  //         max_scores[i][j] = 0.4; //for sigma =1
  //         continue;
  //       }
  //       for(int c = 0; c<=ceil((max-min)*1./stepsize); c++){
  //         float d = min+stepsize*c;
  //         float score_ij = CalculateScoreIJ(d, data,i,j);
  //         if(score_ij > max_score_ij) max_score_ij = score_ij;
  //       }       
  //       max_scores.at(i).at(j) = max_score_ij;
  //     }
  //   }
  //   return max_scores;    
  // }


  template <typename T>
  void WriteVecToFile(std::vector<T>& v, std::ofstream& out_stream){
    unsigned short size = v.size();
    out_stream.write(reinterpret_cast<char*>(&size),sizeof(unsigned short));
    if(size > 0){
      out_stream.write(reinterpret_cast<char*>(&v[0]),size* sizeof(T));
    }
  }


  void SaveDCData(DCData& data, const String& filename){
    std::ofstream out_stream(filename.c_str(), std::ios::binary);
    std::stringstream ss;
    if(! out_stream){
      ss<<"Could not write to file '";
      ss<<filename<<"'";
      throw ost::io::IOException(ss.str());
    }  

    WriteVecToFile(data.cluster_seq_sim, out_stream);
    WriteVecToFile(data.cluster_seq_ids, out_stream);

    unsigned short size = data.tpl_distances.size();
    out_stream.write(reinterpret_cast<char*>(&size),sizeof(unsigned short));
    for(int i = 0; i<size;++i){
      for(int j = i+1; j<size;++j){
        unsigned short size2 = data.get_distances(i,j).size();
        std::vector<unsigned short> short_dist(size2);
        for(int k= 0; k< size2; k++){
          short_dist[k] = data.get_distances(i,j)[k]*100; 
        }
        WriteVecToFile(short_dist, out_stream);
      }
    }

    size = data.cluster_sizes.size();
    out_stream.write(reinterpret_cast<char*>(&size),sizeof(unsigned short));
    for(int i = 0; i<size;++i){
      for(int j = i+1; j<size;++j){
        WriteVecToFile(data.cluster_sizes.at(i).at(j),out_stream);
      }
    }

    size_t ulen = data.seq.length() + 1;
    out_stream.write((const char*)&ulen, sizeof(size_t));
    if(ulen>0) out_stream.write((const char*)data.seq.c_str(), ulen);

    out_stream.write(reinterpret_cast<char*>(&data.dist_cutoff),sizeof(unsigned short));
    out_stream.write(reinterpret_cast<char*>(&data.exp_factor),sizeof(unsigned short));
    out_stream.close();
  }


  template <typename T>
  void LoadVecFromFile(std::vector<T>& v, std::ifstream& in_stream){
    unsigned short size;
    in_stream.read(reinterpret_cast<char*>(&size),sizeof(unsigned short));  
    v.resize(size);
    if(size > 0){
      in_stream.read(reinterpret_cast<char*>(&v[0]),size * sizeof(T));
    }
  }


  DCData LoadDCData(const String& filename){
    std::vector<float> cluster_weights;
    std::vector<float> cluster_seq_sim;
    std::vector<float> cluster_seq_ids;
    std::vector<float> short_dist;
    std::vector<unsigned short> v;
    std::ifstream in_stream(filename.c_str(), std::ios::binary);

    if(! in_stream){
      std::stringstream ss;
      ss<<"Could not read file. File '";
      ss<<filename<<"' does not exist!";
      throw ost::io::IOException(ss.str());
    }
    
    unsigned short size;
    LoadVecFromFile(cluster_seq_sim, in_stream);
    LoadVecFromFile(cluster_seq_ids, in_stream);
      
    in_stream.read(reinterpret_cast<char*>(&size),sizeof(unsigned short));
    std::vector<std::vector<std::vector<float> > > tpl_distances(size);
    for(std::vector<std::vector<std::vector<float> > >::iterator i = tpl_distances.begin(); i != tpl_distances.end(); ++i){
      i->resize(size);
    } 
    for(int i = 0; i<size;++i){
      for(int j = i+1; j<size;++j){
        std::vector<unsigned short> short_dist;
        LoadVecFromFile(short_dist, in_stream);
        std::vector<float> float_dist(short_dist.size());
        for(int k = 0;  k< short_dist.size(); k++){
          float_dist[k] = short_dist[k]/100.;
        }
        tpl_distances.at(i).at(j) = float_dist;
      }
    }
    
    in_stream.read(reinterpret_cast<char*>(&size),sizeof(unsigned short));
    std::vector<std::vector<std::vector<unsigned short> > > cluster_sizes(size);
    for(std::vector<std::vector<std::vector<unsigned short> > >::iterator i = cluster_sizes.begin(); i != cluster_sizes.end(); ++i){
      i->resize(size);
    }
    for(int i = 0; i<size;++i){
      for(int j = i+1; j<size;++j){
        LoadVecFromFile(v, in_stream);
        cluster_sizes.at(i).at(j) = v;
      }
    }

    size_t len = 0;
    in_stream.read((char*)&len, sizeof(size_t));
    char* seq = new char[len];
    if(len>0) in_stream.read(seq, len);

    unsigned short dist_cutoff;
    in_stream.read(reinterpret_cast<char*>(&dist_cutoff),sizeof(unsigned short));  

    unsigned short exp_factor;
    in_stream.read(reinterpret_cast<char*>(&exp_factor),sizeof(unsigned short));  
    in_stream.close();

    for(std::vector<float>::iterator i = cluster_seq_sim.begin(); i != cluster_seq_sim.end(); ++i){
      float weight = exp(exp_factor*(*i));
      cluster_weights.push_back(weight);
    }

    DCData data;
    data.seq = seq;
    data.cluster_weights = cluster_weights;
    data.cluster_seq_sim = cluster_seq_sim;
    data.cluster_seq_ids = cluster_seq_ids; 
    data.tpl_distances = tpl_distances;
    data.cluster_sizes = cluster_sizes;
    data.dist_cutoff = dist_cutoff;
    data.exp_factor = exp_factor;
    return data;
  }


  DCData FillDCData(const ost::seq::AlignmentHandle& msaln, const std::vector<std::vector<int> >& cluster, const std::vector<String>& filenames,
                    ost::seq::alg::SubstWeightMatrixPtr subst, unsigned short dist_cutoff /*=15*/, unsigned short exp_factor /*=70*/){

    if(msaln.GetCount()-1 != filenames.size()){
      throw std::runtime_error("Length of filename list must match the number of templates in msaln");
    }

    if(msaln.GetCount() == 0){
      throw std::runtime_error("msaln must contain a sequence");
    }    

    std::vector<std::vector<geom::Vec3> > structural_info;
    std::vector<bool> valid_structures;
    std::vector<float> cluster_weights;
    std::vector<float> cluster_seq_sim;
    std::vector<float> cluster_seq_ids;
    std::vector<std::vector<std::vector<float> > >  tpl_distances;
    std::vector<std::vector<std::vector<unsigned short> > > cluster_sizes;

    structural_info = ReadCalphaPos(filenames, valid_structures); 
    cluster_seq_sim = CalculateClusterSeqSim(msaln, cluster, subst, valid_structures);
    for(std::vector<float>::iterator i = cluster_seq_sim.begin(); i != cluster_seq_sim.end(); ++i){
      float weight = exp(exp_factor*(*i));
      cluster_weights.push_back(weight);
    }
    
    tpl_distances = ExtractDistances(msaln,cluster,structural_info, valid_structures, cluster_sizes, dist_cutoff);
    cluster_seq_ids = CalculateClusterSeqIds(msaln,cluster, valid_structures);

    DCData data;
    data.seq = msaln.GetSequence(0).GetGaplessString();
    data.cluster_sizes = cluster_sizes;
    data.tpl_distances = tpl_distances;
    data.cluster_seq_sim = cluster_seq_sim;
    data.cluster_weights = cluster_weights;
    data.cluster_seq_ids = cluster_seq_ids;
    data.dist_cutoff = dist_cutoff;
    data.exp_factor = exp_factor;
    return data;  
  }


  std::vector<std::vector<float> > DetermineFeatureValues(const ost::seq::AlignmentHandle& aln, const DCData& data){

    int cluster_start, cluster_end, j, i=0;
    float dist;
    const unsigned short dist_cutoff = data.dist_cutoff;
    const int num_features = 5;
    const int aln_length = aln.GetLength();
    float feature_values[aln_length][num_features];   
    memset(feature_values, 0, sizeof(feature_values[0][0]) * aln_length * num_features);     
    std::vector<std::vector<float> > model_feature_values; 

    //check for valid DCData and valid alignment
    if(data.tpl_distances.size() != aln.GetLength())  throw std::invalid_argument("DCData does not match alignment!");
    if(aln.GetSequence(0).GetString()!= data.seq.c_str() ){
      throw std::invalid_argument("SEQRES of DCData does not match SEQRES in alignment!");
    }

    //Compute features
    for (ost::seq::AlignmentHandle::iterator aln_it1 = aln.begin(), e1 = aln.end(); aln_it1 != e1; aln_it1++, i++) {
      const ost::seq::AlignedColumn& col1 = *aln_it1;
      j = i+1;
      if (col1[1] == '-') continue;
      for (ost::seq::AlignmentHandle::iterator aln_it2 = ost::seq::AlignedColumnIterator(aln,j,aln.GetLength()), e2 = aln.end(); aln_it2 != e2; aln_it2++, j++) {            
        const ost::seq::AlignedColumn& col2 = *aln_it2;
        if (col2[1] == '-') continue;
        ost::mol::ResidueView res_i = aln.GetSequence(1).GetResidue(i);
        ost::mol::AtomView ca_i = res_i.FindAtom("CA");
        ost::mol::ResidueView res_j = aln.GetSequence(1).GetResidue(j);             
        ost::mol::AtomView ca_j = res_j.FindAtom("CA"); 
        if(!ca_i.IsValid() || !ca_j.IsValid()) continue; 
        dist = geom::Distance(ca_i.GetPos(),ca_j.GetPos());
        if(dist >= dist_cutoff) continue;
        if(data.get_distances(i,j).empty()) continue;
        std::vector<float> distances = data.get_distances(i,j);
        std::vector<unsigned short> cluster_idx = data.get_cluster(i,j);

        float max_seq_id  = 0;
        float num_cluster = 0;
        float max_seq_sim = 0;
        cluster_start = 0;
        for(int k = 0; k< cluster_idx.size(); k++){          
          cluster_end = cluster_idx[k];
          if(cluster_start == cluster_end) continue;
          if(data.cluster_seq_ids.at(k)> max_seq_id) max_seq_id = data.cluster_seq_ids[k];
          if(data.cluster_seq_sim[k]> max_seq_sim) max_seq_sim = data.cluster_seq_sim[k];
          num_cluster ++;
          cluster_start = cluster_end;
        }

        float var = CalculateVariance(distances);
        int num_tpls = distances.size();

        feature_values[i][0] += max_seq_sim;
        feature_values[j][0] += max_seq_sim;
        feature_values[i][1] += num_cluster;
        feature_values[j][1] += num_cluster;
        feature_values[i][2] += var;
        feature_values[j][2] += var;
        feature_values[i][3] += 1;  //count residues within cutoff radius
        feature_values[j][3] += 1;
        feature_values[i][4] += max_seq_id;
        feature_values[j][4] += max_seq_id;        
        // feature_values[i][5] += num_tpls;
        // feature_values[j][5] += num_tpls;
      }
    }

    int first_res = aln.GetSequence(1).GetFirstNonGap();
    int last_res = aln.GetSequence(1).GetLastNonGap();   
    for(int i = first_res; i<=last_res; i++ ){
      std::vector<float> res_feature_values;
      if(aln.GetSequence(1)[i] == '-') continue;
      for(int j = 0; j< num_features; j++){
        float num_inter = feature_values[i][3];
        if(num_inter == 0)  res_feature_values.push_back(NAN);
        else{
          if(j == 3) res_feature_values.push_back(num_inter);
          else res_feature_values.push_back(feature_values[i][j]/(num_inter));
        }
      } 
      model_feature_values.push_back(res_feature_values);
    }
    return model_feature_values;
  }


  float CalculateVariance(const std::vector<float>& distances){
  
    float size = distances.size();
    float mean = std::accumulate(distances.begin(), distances.end(), 0.0)/size;
    float temp = 0;
    for(int i = 0; i < size; i++){
      temp += (distances[i] - mean)*(distances[i] - mean);
    }
    return temp/size;
  }


}

