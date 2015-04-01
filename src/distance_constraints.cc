#include <qmean/distance_constraints.hh>
#include <ost/seq/aligned_column.hh>
#include <ost/seq/aligned_column_iterator.hh>




namespace qmean{

std::map<IndexPair, float >  DistanceConstraints(const ost::seq::AlignmentHandle& aln) {
  std::map<IndexPair, float > mus;

  const int s = 2;  
  const float D = 15; 
  float mu_ijk;
  int idx0, idx1, k2, k1 = 0;
  IndexPair p;

  for (ost::seq::AlignmentHandle::iterator i = aln.begin(), e1 = aln.end(); i != e1; i++, k1++) {
    const ost::seq::AlignedColumn& col1 = *i;
    k2 = k1+s;
    if ((col1[0] != '-') && (col1[1] != '-') && k1+s < aln.GetLength()){
       
      for (ost::seq::AlignmentHandle::iterator j = ost::seq::AlignedColumnIterator(aln,k2,aln.GetLength()), e = aln.end(); j != e; j++, k2++) {
          
        const ost::seq::AlignedColumn& col2 = *j;
        if ((col2[0] != '-') && (col2[1] != '-')){
          ost::mol::ResidueView res1 = aln.GetSequence(1).GetResidue(k1);
          ost::mol::AtomView ca1 = res1.FindAtom("CA");
          ost::mol::ResidueView res2 = aln.GetSequence(1).GetResidue(k2);             
          ost::mol::AtomView ca2 = res2.FindAtom("CA");           
          if(ca1.IsValid() && ca2.IsValid()){            
            mu_ijk = geom::Distance(ca1.GetPos(),ca2.GetPos());


            if(mu_ijk < D){
              idx0 = aln.GetResidueIndex(0,k1);
              idx1 = aln.GetResidueIndex(0,k2);

              p.index_one = idx0;
              p.index_two = idx1; 

              mus[p] = mu_ijk;
            }
          }
        }
      }
    }
  }
  return  mus;
}









float seq_sim(const ost::seq::alg::SubstWeightMatrix& mat, const ost::seq::ConstSequenceHandle& seq1, const ost::seq::ConstSequenceHandle& seq2){
  float seq_sim, sim_score= 0;
  int length= 0;
  const int min_s = -4;
  const int max_s = 11;

  int k = 0;
  for (int i = 0; i<seq1.GetLength(); ++i) {  
    if( seq1[i]!= '-' && seq2[i]!='-'){
      sim_score+= mat.GetWeight(seq1[i], seq2[i]);
      length+= 1;
    }
  }    
  if(length != 0){
    seq_sim = sim_score/float(length); 
    seq_sim = (seq_sim-min_s)/float(max_s-min_s);
  }
  else{
    seq_sim = 0;
  }
  return seq_sim;
}





}



















