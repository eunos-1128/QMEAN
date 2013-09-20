#include <qmean/spherical_smoother.hh>


namespace qmean{

SphericalSmoother::SphericalSmoother(std::vector<geom::Vec3>& positions, Real std){

  ost::mol::SpatialOrganizer<SpatialOrganizerItemSmoother, geom::Vec3>  spatial_organizer(5);
  
  //build the spatial organizer and build up the value distributor
  std::vector<geom::Vec3>::iterator pos_it=positions.begin();
  int index=0;
  for(; pos_it!=positions.end();++pos_it,++index){
    spatial_organizer.Add(SpatialOrganizerItemSmoother(*pos_it,index),*pos_it);
  }  

  pos_it=positions.begin();
  Real weight;
  Real dist;
  Real squared_std=std*std;
  index=0;

  for(; pos_it!=positions.end();++pos_it,++index){

    std::vector<std::pair<Real,Real> > visits;
    std::vector<int> visit_indices;
    std::vector<SpatialOrganizerItemSmoother> in_reach = spatial_organizer.FindWithin(*pos_it,3*std);

    std::vector<SpatialOrganizerItemSmoother>::iterator reach_it=in_reach.begin();
    for(;reach_it!=in_reach.end();++reach_it){
      dist=geom::Distance(reach_it->pos,*pos_it);
      weight=exp(-0.5*(dist*dist/squared_std));
      std::pair<Real,Real> weight_value_pair=std::make_pair(weight,0.0);
      visits.push_back(weight_value_pair);
      visit_indices.push_back(reach_it->index);
    }
    visitor_pattern_.push_back(visits);
    value_distributor_.push_back(visit_indices);
  }

}

std::vector<Real> SphericalSmoother::Smooth(std::vector<Real>& values){

  if(values.size()!=visitor_pattern_.size()){
    std::stringstream ss;
    ss << "Value list has not the same length as the position list from which the smoother object has been created!";
    throw ost::io::IOException(ss.str());
  }

  //fill the values in the visitor pattern
  int j;
  for(int i=0;i<values.size();++i){
    j=0;
    for(std::vector<int>::iterator it=value_distributor_[i].begin();it!=value_distributor_[i].end();++it,++j){
      visitor_pattern_[i][j].second=values[*it];
    }
  }

  //create and return the result

  std::vector<Real> result;

  for(int i=0;i<values.size();++i){
    //check for nan
    if(values[i]!=values[i]){
      result.push_back(std::numeric_limits<Real>::quiet_NaN());
      continue;
    }
    result.push_back(this->Evaluate(visitor_pattern_[i]));
  }

  return result;
}

Real SphericalSmoother::Evaluate(std::vector<std::pair<Real,Real> >& pair_vec){
  std::vector<Real> weights;
  std::vector<Real> values;
  std::vector<std::pair<Real,Real> >::iterator it=pair_vec.begin();
  for(;it!=pair_vec.end();++it){
    //check for nan
    if(it->second!=it->second or it->first!=it->first){
      continue;
    }
    weights.push_back(it->first);
    values.push_back(it->second);
  }

  if(values.size()==0){
    //no valid values at all!
    return std::numeric_limits<Real>::quiet_NaN();
  }

  //normalize weights
  Real sum_of_weights=0.0;
  for(std::vector<Real>::iterator it=weights.begin(); it!=weights.end(); ++it){
    sum_of_weights+=(*it);
  }

  for(std::vector<Real>::iterator it=weights.begin(); it!=weights.end(); ++it){
    (*it)/=sum_of_weights;
  }

  //calculate final value
  Real final_value = 0.0;
  std::vector<Real>::iterator w_it=weights.begin();
  std::vector<Real>::iterator v_it=values.begin();
  
  for(;w_it!=weights.end(); ++w_it, ++v_it){
    final_value+=((*w_it)*(*v_it));
  }
  return final_value;
}

}//namespace
