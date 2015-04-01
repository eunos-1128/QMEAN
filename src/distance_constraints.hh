#include <ost/geom/geom.hh>
#include <ost/mol/mol.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/subst_weight_matrix.hh>




namespace qmean{

struct IndexPair{
  IndexPair(): index_one(0),index_two(0) { }
  IndexPair(int one,int two): index_one(one), index_two(two) { }
  int index_one;
  int index_two;


  bool operator < (const IndexPair& other) const{
    return ( index_one < other.index_one ) ||
           (( index_one == other.index_one) && 
            ( index_two < other.index_two ));
  }
};

inline std::ostream&
operator<<(std::ostream& o, const IndexPair& g){
   o << '(' << g.index_one << ',' << g.index_two << ')';
   return o;
}





std::map<IndexPair, float > DistanceConstraints(const ost::seq::AlignmentHandle& aln);
float seq_sim(const ost::seq::alg::SubstWeightMatrix& mat, const ost::seq::ConstSequenceHandle& seq1, const ost::seq::ConstSequenceHandle& seq2);





}
