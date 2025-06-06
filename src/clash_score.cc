// Copyright (c) 2013-2020, SIB - Swiss Institute of Bioinformatics and
// Biozentrum - University of Basel
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <ost/mol/atom_view.hh>
#include <ost/mol/residue_view.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/chain_view.hh>
#include <ost/mol/entity_view.hh>
#include <ost/mol/bond_handle.hh>
#include <ost/geom/vecmat3_op.hh>
#include <qmean/clash_score.hh>

namespace qmean {

std::vector<Real> GetClashScores(const ost::mol::EntityView& target, 
                                 const ost::mol::EntityView& env) {

  if(target.GetHandle() != env.GetHandle()) {
    throw std::runtime_error("target and env are expected to be of the same "
                             "handle");
  }

  ost::mol::ResidueViewList res_list = target.GetResidueList();
  std::vector<Real> return_scores(res_list.size(), 0.0);
  
  for(uint r_idx = 0; r_idx < res_list.size(); ++r_idx) {
    
    ost::mol::AtomViewList atom_list = res_list[r_idx].GetAtomList();
    int rnum = res_list[r_idx].GetNumber().GetNum(); 
    String chain_name = res_list[r_idx].GetChain().GetName();
    
    for(ost::mol::AtomViewList::const_iterator at_it = atom_list.begin();
        at_it != atom_list.end(); ++at_it) {
 
      String aname = at_it->GetName();
      if(aname.empty()) {
        continue;
      }

      ost::mol::AtomViewList close_atoms = env.FindWithin(at_it->GetPos(), 
                                                          Real(3.2));
      
      for(ost::mol::AtomViewList::const_iterator close_at_it = close_atoms.begin();
          close_at_it != close_atoms.end(); ++close_at_it) {

        // we're not assessing interactions from the same or neighboring residues
        int close_rnum = close_at_it->GetResidue().GetNumber().GetNum();
        if(std::abs(close_rnum - rnum) < 2){
          if(chain_name == close_at_it->GetResidue().GetChain().GetName()){
            continue;
          }
        }

        String close_aname = close_at_it->GetName();
        if(close_aname.empty()) {
          continue;
        }

        // we're  not assessing disulfid bridges this also includes interactions
        // of gamma sulfurs towards cbeta atoms of other cysteins
        if(aname == "SG") {
          if(close_aname == "SG" || 
             (close_at_it->GetResidue().GetName() == "CYS" && 
              close_at_it->GetName() == "CB")) {
            continue;
          }
        }
        if(close_aname == "SG") {
          if(at_it->GetResidue().GetName() == "CYS" && at_it->GetName() == "CB") {
            continue;
          }
        }


        Real d = geom::Distance(at_it->GetPos(), close_at_it->GetPos());

        // fetching the element by atom name obviously only works if we're dealing
        // with peptide only views... But thats an assumption we're doing in QMEAN
        // anyway...
        char ele1 = aname[0];
        char ele2 = close_aname[0];
        Real r1 = 0.0;
        Real r2 = 0.0;

        if(ele1 == 'C') {
          r1 = 1.6;
        }
        else if(ele1 == 'N') {
          r1 = 1.3;
        }
        else if(ele1 == 'O') {
          r1 = 1.3;
        } 
        else if(ele1 == 'S') {
          r1 = 1.7;
        }

        if(ele2 == 'C') {
          r2 = 1.6;
        }
        else if(ele2 == 'N') {
          r2 = 1.3;
        }
        else if(ele2 == 'O') {
          r2 = 1.3;
        } 
        else if(ele2 == 'S') {
          r2 = 1.7;
        }

        Real clash_score = 0.0;

        if(r1 != 0.0 && r2 != 0.0) {
          Real Rij = r1 + r2;
          if(d < Real(0.8254) * Rij) {
            clash_score = 10.0;
          } else if(d <= Rij) {
            clash_score = 57.273 * (1 - d/Rij);
          }
        }

        return_scores[r_idx] += clash_score;
      }
    }
  }

  return return_scores;
}

} // ns
