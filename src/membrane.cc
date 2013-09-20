#include <qmean/membrane.hh>
#include <ctime>
#include <cstdlib>
#include <sys/time.h>

SolvationGrid::SolvationGrid(geom::AlignedCuboid cub, Real bin_size):
  bin_size_(bin_size){

  extent_ = cub.GetMax()-cub.GetMin();

  xbins_ = int(extent_[0]/bin_size);
  ybins_ = int(extent_[1]/bin_size);
  zbins_ = int(extent_[2]/bin_size);

  origin_ = cub.GetMin();

  for(int i=0;i<xbins_;++i){
    std::vector<std::vector<int> > two_d;
    for(int j=0;j<ybins_;++j){
      std::vector<int> one_d;
      for(int k=0;k<zbins_;++k){
        one_d.push_back(0);
      }
      two_d.push_back(one_d);
    }
    grid_.push_back(two_d); 
  }
  //create invalid views
  surf_ = ost::mol::SurfaceHandle();
  solvated_surf_ = ost::mol::SurfaceHandle();

  is_flooded_=false;

}

ost::img::ImageHandle SolvationGrid::AsIMG(){

  ost::img::ImageHandle img = ost::img::CreateImage(ost::img::Size(xbins_,ybins_,zbins_));
  img.SetAbsoluteOrigin(origin_);
  img.SetPixelSampling(bin_size_);

  for(int i=0;i<xbins_;++i){
    for(int j=0;j<ybins_;++j){
      for(int k=0;k<zbins_;++k){
        img.SetReal(ost::img::Point(i,j,k),grid_[i][j][k]);
      }
    }
  }

  return img;
}

void SolvationGrid::AddSurface(const ost::mol::SurfaceHandle& surf){

  surf_ = surf;
  ost::mol::SurfaceTriIDList tri_id = surf.GetTriIDList();
  ost::mol::SurfaceTri tri;

  int i,j,k;

  geom::Vec3 u0;
  geom::Vec3 u1;
  geom::Vec3 u2;
  geom::Vec3 cube_center;

  Real min_x;
  Real max_x;
  Real min_y;
  Real max_y;
  Real min_z;
  Real max_z;

  int min_x_idx;
  int max_x_idx;
  int min_y_idx;
  int max_y_idx;
  int min_z_idx;
  int max_z_idx;

  Real half_bin_size = bin_size_/2;

  for(ost::mol::SurfaceTriIDList::iterator it = tri_id.begin(); it!=tri_id.end(); ++it){

    tri = surf.GetTri(*it);
    u0 = surf.GetVertex(tri.v0).position;
    u1 = surf.GetVertex(tri.v1).position;
    u2 = surf.GetVertex(tri.v2).position;
    //extract region around triangle
    min_x = std::min(u0[0],std::min(u1[0],u2[0]));
    min_y = std::min(u0[1],std::min(u1[1],u2[1]));
    min_z = std::min(u0[2],std::min(u1[2],u2[2]));

    max_x = std::max(u0[0],std::max(u1[0],u2[0]));
    max_y = std::max(u0[1],std::max(u1[1],u2[1]));
    max_z = std::max(u0[2],std::max(u1[2],u2[2]));

    min_x_idx = std::floor((min_x-origin_[0])/bin_size_);
    max_x_idx = std::ceil((max_x-origin_[0])/bin_size_);
    min_y_idx = std::floor((min_y-origin_[1])/bin_size_);
    max_y_idx = std::ceil((max_y-origin_[1])/bin_size_);
    min_z_idx = std::floor((min_z-origin_[2])/bin_size_);
    max_z_idx = std::ceil((max_z-origin_[2])/bin_size_);

    for(i=std::max(min_x_idx,0);i<=std::min(max_x_idx,xbins_-1);++i){
      for(j=std::max(min_y_idx,0);j<=std::min(max_y_idx,ybins_-1);++j){
        for(k=std::max(min_z_idx,0);k<=std::min(max_z_idx,zbins_-1);++k){
          cube_center[0]=origin_[0]+i*bin_size_+half_bin_size;
          cube_center[1]=origin_[1]+j*bin_size_+half_bin_size;
          cube_center[2]=origin_[2]+k*bin_size_+half_bin_size;
          if(this->InCube(cube_center,u0,u1,u2)){
            grid_[i][j][k] = 1;
          }
        }
      }
    }
  }

}

bool SolvationGrid::InCube(const geom::Vec3& cube_center,const geom::Vec3& u0,const geom::Vec3& u1, const geom::Vec3& u2){
  

  geom::Vec3 v0 = u0-cube_center;
  geom::Vec3 v1 = u1-cube_center;
  geom::Vec3 v2 = u2-cube_center;

  Real p0;
  Real p1;
  Real p2;
  Real half_bin_size = bin_size_/2;
  Real r;

  //check normals on cube (three cases)

  //e0
  p0 = v0[0];
  p1 = v1[0];
  p2 = v2[0];
  r = half_bin_size;

  if(std::min(p0,std::min(p1,p2))>r || std::max(p0,std::max(p1,p2))<-r){
    return false;
  }

  //e1
  p0 = v0[1];
  p1 = v1[1];
  p2 = v2[1];

  if(std::min(p0,std::min(p1,p2))>r || std::max(p0,std::max(p1,p2))<-r){
    return false;
  }
  
  //e2
  p0 = v0[2];
  p1 = v1[2];
  p2 = v2[2];

  if(std::min(p0,std::min(p1,p2))>r || std::max(p0,std::max(p1,p2))<-r){
    return false;
  }

  //check cross products of edges (9 cases)

  //a00
  p0 = v0[2]*v1[1]-v0[1]*v1[2];
  p2 = v2[2]*(v1[1]-v0[1])+v2[1]*(v0[2]-v1[2]);
  r = half_bin_size*(std::abs(v0[2]-v1[2])+std::abs(v1[1]-v0[1]));
  if(std::min(p0,p2)>r || std::max(p0,p2)<-r){
    return false;
  }

  //a01
  p0 = v0[1]*(v1[2]-v2[2])+v0[2]*(v2[1]-v1[1]);
  p1 = v1[2]*v2[1]-v1[1]*v2[2];
  r = half_bin_size*(std::abs(v1[2]-v2[2])+std::abs(v2[1]-v1[1]));
  if(std::min(p0,p1)>r || std::max(p0,p1)<-r){
    return false;
  }

  //a02
  p0 = v0[1]*v2[2]-v0[2]*v2[1];
  p1 = v1[1]*(v2[2]-v0[2])+v1[2]*(v0[1]-v2[1]);
  r = half_bin_size*(std::abs(v2[2]-v0[2])+std::abs(v0[1]-v2[1]));
  if(std::min(p0,p1)>r || std::max(p0,p1)<-r){
    return false;
  }

  //a10
  p0 = v0[0]*v1[2]-v0[2]*v1[0];
  p2 = v2[2]*(v0[0]-v1[0])+v2[0]*(v1[2]-v0[2]);
  r = half_bin_size*(std::abs(v1[2]-v0[2])+std::abs(v0[0]-v1[0]));
  if(std::min(p0,p2)>r || std::max(p0,p2)<-r){
    return false;
  }

  //a11
  p0 = v0[0]*(v2[2]-v1[2])+v0[2]*(v1[0]-v2[0]);
  p1 = v1[0]*v2[2]-v1[2]*v2[0];
  r = half_bin_size*(std::abs(v2[2]-v1[2])+std::abs(v1[0]-v2[0]));
  if(std::min(p0,p1)>r || std::max(p0,p1)<-r){
    return false;
  }

  //a12
  p0 = v0[2]*v2[0]-v0[0]*v2[2];
  p1 = v1[0]*(v0[2]-v2[2])+v1[2]*(v2[0]-v0[0]);
  r = half_bin_size*(std::abs(v0[2]-v2[2])+std::abs(v2[0]-v0[0]));
  if(std::min(p0,p1)>r || std::max(p0,p1)<-r){
    return false;
  }

  //a20
  p0 = v0[1]*v1[0]-v0[0]*v1[1];
  p2 = v2[1]*(v1[0]-v0[0])+v2[0]*(v0[1]-v1[1]);
  r = half_bin_size*(std::abs(v0[1]-v1[1])+std::abs(v1[0]-v0[0]));
  if(std::min(p0,p2)>r || std::max(p0,p2)<-r){
    return false;
  }

  //a21
  p0 = v0[0]*(v1[1]-v2[1])+v0[1]*(v2[0]-v1[0]);
  p1 = v1[1]*v2[0]-v1[0]*v2[1];
  r = half_bin_size*(std::abs(v1[1]-v2[1])+std::abs(v2[0]-v1[0]));
  if(std::min(p0,p1)>r || std::max(p0,p1)<-r){
    return false;
  }

  //a22
  p0 = v0[0]*v2[1]-v0[1]*v2[0];
  p1 = v1[0]*(v2[1]-v0[1])+v1[1]*(v0[0]-v2[0]);
  r = half_bin_size*(std::abs(v2[1]-v0[1])+std::abs(v0[0]-v2[0]));
  if(std::min(p0,p1)>r || std::max(p0,p1)<-r){
    return false;
  }

  //check normal of triangle

  p0 = v2[2]*(v0[0]*v1[1]-v0[1]*v1[0])+v2[1]*(v0[2]*v1[0]-v0[0]*v1[2])+v2[0]*(v0[1]*v1[2]-v0[2]*v1[1]);
  r = half_bin_size*(std::abs((v1[1]-v0[1])*(v2[2]-v0[2])-(v1[2]-v0[2])*(v2[1]-v0[1])));
  r += half_bin_size*(std::abs((v1[2]-v0[2])*(v2[0]-v0[0])-(v1[0]-v0[0])*(v2[2]-v0[2])));
  r += half_bin_size*(std::abs((v1[0]-v0[0])*(v2[1]-v0[1])-(v1[1]-v0[1])*(v2[0]-v0[0])));
  if(p0 > r || p0 < -r){
    return false;
  }

  return true;

}



void SolvationGrid::Flood(){
  for(int level=0;level<zbins_;++level){
    this->FloodLevel(level,std::make_pair(0,0),0,2);
    this->FloodLevel(level,std::make_pair(xbins_-1,ybins_-1),0,2);
    this->FloodLevel(level,std::make_pair(xbins_-1,0),0,2);
    this->FloodLevel(level,std::make_pair(0,ybins_-1),0,2);
  }
  is_flooded_=true;
}



void SolvationGrid::FloodLevel8Neighbours(int level, std::pair<int,int> start_coordinates, int orig_value, int dest_value){
 
  if(orig_value!=grid_[start_coordinates.first][start_coordinates.second][level]){
    return;
  }

  std::vector<std::pair<int,int> > queue;

  queue.push_back(std::make_pair(start_coordinates.first,start_coordinates.second));

  std::pair<int,int> actual_position;

  while(!queue.empty()){
    actual_position = queue.back();
    queue.pop_back();
    //check xrange
    if(actual_position.first<0 || actual_position.first>=xbins_){
      continue;
    }
    //check yrange
    if(actual_position.second<0 || actual_position.second>=ybins_){
      continue;
    }
    //check whether there is an atom
    if(grid_[actual_position.first][actual_position.second][level]!=orig_value){
      continue;
    }
    //check whether we already set the value
    if(grid_[actual_position.first][actual_position.second][level]==dest_value){
      continue;
    }

    //set value
    grid_[actual_position.first][actual_position.second][level]=dest_value;
    

    queue.push_back(std::make_pair(actual_position.first,actual_position.second+1));
    queue.push_back(std::make_pair(actual_position.first,actual_position.second-1));
    queue.push_back(std::make_pair(actual_position.first+1,actual_position.second));
    queue.push_back(std::make_pair(actual_position.first-1,actual_position.second));
    queue.push_back(std::make_pair(actual_position.first-1,actual_position.second-1));
    queue.push_back(std::make_pair(actual_position.first+1,actual_position.second+1));
    queue.push_back(std::make_pair(actual_position.first+1,actual_position.second-1));
    queue.push_back(std::make_pair(actual_position.first-1,actual_position.second+1));
  }
}


void SolvationGrid::FloodLevel(int level, std::pair<int, int> start_coordinates, int orig_value, int dest_value){
  //http://lodev.org/cgtutor/floodfill.html

  if(orig_value!=grid_[start_coordinates.first][start_coordinates.second][level]){
    return;
  }

  std::vector<std::pair<int,int> > queue;
  queue.push_back(start_coordinates);


    
  int y1,y,x; 
  bool spanLeft, spanRight;
  std::pair<int,int> actual_position;
    
  while(!queue.empty()){    

    actual_position = queue.back();
    queue.pop_back();
    x = actual_position.first;
    y = actual_position.second;
    y1 = y;

    while(y1 >= 0 && grid_[x][y1][level] == orig_value) y1--;
    y1++;
    spanLeft = spanRight = 0;
    while(y1 < grid_[0].size() && grid_[x][y1][level] == orig_value )
    {
      grid_[x][y1][level] = dest_value;
      if(!spanLeft && x > 0 && grid_[x-1][y1][level] == orig_value) 
      {
        queue.push_back(std::make_pair(x-1,y1));
        spanLeft = 1;
      }
      else if(spanLeft && x > 0 && grid_[x-1][y1][level] != orig_value)
      {
        spanLeft = 0;
      }
      if(!spanRight && x < grid_.size() - 1 && grid_[x+1][y1][level] == orig_value) 
      {
        queue.push_back(std::make_pair(x+1,y1));
        spanRight = 1;
      }
      else if(spanRight && x < grid_.size() - 1 && grid_[x+1][y1][level] != orig_value)
      {
        spanRight = 0;
      } 
      y1++;
    }
  }
}

bool SolvationGrid::IsFilled(geom::Vec3& pos){

  int x_index = int(std::floor((pos[0]-origin_[0])/bin_size_));
  int y_index = int(std::floor((pos[1]-origin_[1])/bin_size_));
  int z_index = int(std::floor((pos[2]-origin_[2])/bin_size_)); 

  if(x_index<0 || x_index>=xbins_ || y_index<0 || y_index>=ybins_ || z_index<0 || z_index>=zbins_){
    return false;
  }

  return grid_[x_index][y_index][z_index]==2;
}


void SolvationGrid::CalculateSolvatedSurface(){
  if(!surf_.IsValid()){
    throw "You have to add surface first!";
  }

  if(!is_flooded_){
    this->Flood();
    is_flooded_=true;
  }

  solvated_surf_ = ost::mol::CreateSurface();

  int level,x,y;

  for(level=0;level<zbins_;++level){
    for(x=0;x<xbins_;++x){
      for(y=0;y<ybins_;++y){
        if(grid_[x][y][level]!=2){
          this->FloodLevel(level, std::make_pair(x,y), 1, 10);
          break;
        }
      }
      for(y=ybins_-1;y>=0;--y){
        if(grid_[x][y][level]!=2){
          this->FloodLevel(level, std::make_pair(x,y), 1, 10);
          break;
        }
      }
    }
    for(y=0;y<ybins_;++y){
      for(x=0;x<xbins_;++x){
        if(grid_[x][y][level]!=2){
          this->FloodLevel(level, std::make_pair(x,y), 1, 10);
          break;
        }
      }
      for(x=xbins_-1;x>=0;--x){
        if(grid_[x][y][level]!=2){
          this->FloodLevel(level, std::make_pair(x,y), 1, 10);
          break;
        }
      }
    }
  }


  ost::mol::SurfaceTriIDList tri_id = surf_.GetTriIDList();
  ost::mol::SurfaceTri tri;

  int i,j,k;

  geom::Vec3 u0;
  geom::Vec3 u1;
  geom::Vec3 u2;
  geom::Vec3 cube_center;

  Real min_x;
  Real max_x;
  Real min_y;
  Real max_y;
  Real min_z;
  Real max_z;

  int min_x_idx;
  int max_x_idx;
  int min_y_idx;
  int max_y_idx;
  int min_z_idx;
  int max_z_idx;

  ost::mol::SurfaceVertex v0;
  ost::mol::SurfaceVertex v1;
  ost::mol::SurfaceVertex v2;

  ost::mol::SurfaceVertexID id0;
  ost::mol::SurfaceVertexID id1;
  ost::mol::SurfaceVertexID id2;

  Real half_bin_size = bin_size_/2;

  bool added_vertex;

  for(ost::mol::SurfaceTriIDList::iterator it = tri_id.begin(); it!=tri_id.end(); ++it){

    added_vertex = false;

    tri = surf_.GetTri(*it);
    u0 = surf_.GetVertex(tri.v0).position;
    u1 = surf_.GetVertex(tri.v1).position;
    u2 = surf_.GetVertex(tri.v2).position;
    //extract region around triangle
    min_x = std::min(u0[0],std::min(u1[0],u2[0]));
    min_y = std::min(u0[1],std::min(u1[1],u2[1]));
    min_z = std::min(u0[2],std::min(u1[2],u2[2]));

    max_x = std::max(u0[0],std::max(u1[0],u2[0]));
    max_y = std::max(u0[1],std::max(u1[1],u2[1]));
    max_z = std::max(u0[2],std::max(u1[2],u2[2]));

    min_x_idx = std::floor((min_x-origin_[0])/bin_size_);
    max_x_idx = std::ceil((max_x-origin_[0])/bin_size_);
    min_y_idx = std::floor((min_y-origin_[1])/bin_size_);
    max_y_idx = std::ceil((max_y-origin_[1])/bin_size_);
    min_z_idx = std::floor((min_z-origin_[2])/bin_size_);
    max_z_idx = std::ceil((max_z-origin_[2])/bin_size_);

    for(i=std::max(min_x_idx,0);i<=std::min(max_x_idx,xbins_-1);++i){
      for(j=std::max(min_y_idx,0);j<=std::min(max_y_idx,ybins_-1);++j){
        for(k=std::max(min_z_idx,0);k<=std::min(max_z_idx,zbins_-1);++k){
          cube_center[0]=origin_[0]+i*bin_size_+half_bin_size;
          cube_center[1]=origin_[1]+j*bin_size_+half_bin_size;
          cube_center[2]=origin_[2]+k*bin_size_+half_bin_size;
          if(this->InCube(cube_center,u0,u1,u2)){
            if(grid_[i][j][k] == 10){
              v0 = surf_.GetVertex(tri.v0);
              v1 = surf_.GetVertex(tri.v1);
              v2 = surf_.GetVertex(tri.v2);
              id0 = solvated_surf_.AddVertex(v0);
              id1 = solvated_surf_.AddVertex(v1);
              id2 = solvated_surf_.AddVertex(v2);
              solvated_surf_.AddTri(id0, id1, id2);
              added_vertex = true;
              break;
            }
          }
        }
      
      if(added_vertex){
        break;
      }
    }
    if(added_vertex){
      break;
    }
    }
  }

  for(level=0;level<zbins_;++level){
    for(x=0;x<xbins_;++x){
      for(y=0;y<ybins_;++y){
        if (grid_[x][y][level]==10){
          grid_[x][y][level]=1;
        }
      }
    }
  }
}

ost::mol::SurfaceHandle SolvationGrid::GetSolvatedSurface(){

  if(!solvated_surf_.IsValid()){
    this->CalculateSolvatedSurface();
  }

  return solvated_surf_;

}


void SolvationGrid::CalculateSolvatedView(){
  if(!surf_.IsValid()){
    throw "You have to add surface first!";
  }
  if(!view_.IsValid()){
    throw "You have to add view first!";
  }

  if(!solvated_surf_.IsValid()){
    this->CalculateSolvatedSurface();
  }


  ost::mol::SurfaceVertexIDList v_list = solvated_surf_.GetVertexIDList();
  ost::mol::SurfaceVertex vert;
  ost::mol::AtomViewList a_list;

  for(ost::mol::SurfaceVertexIDList::iterator it=v_list.begin();it!=v_list.end();++it){
    vert = solvated_surf_.GetVertex(*it);
    a_list = view_.FindWithin(vert.position,3.0);
    if(a_list.empty()) continue;
    ost::mol::AtomViewList::const_iterator ait=a_list.begin();
    Real best_dist2 = geom::Length2(ait->GetPos()-vert.position);
    ost::mol::AtomView best_atom = *ait;
    ++ait;
    for(;ait!=a_list.end();++ait) {
      Real dist2 = geom::Length2(ait->GetPos()-vert.position);
      if(dist2<best_dist2) {
        best_dist2=dist2;
        best_atom = *ait;
      }
    }

    best_atom.SetIntProp("surface_atom",1);
  }

  solvated_view_ = view_.Select("gasurface_atom:0=1");
}

ost::mol::EntityView SolvationGrid::GetSolvatedView(){
  if(!solvated_view_.IsValid()){
    this->CalculateSolvatedView();
  }

  return solvated_view_;

}


ost::mol::EntityHandle FillMembraneDummies(const geom::AlignedCuboid& cuboid, const ost::mol::SurfaceHandle& surf, Real solvation_grid_bin_size, Real density){


  SolvationGrid grid(cuboid, solvation_grid_bin_size);

  //grid.AddView(return_view);
  grid.AddSurface(surf);
  grid.Flood(); 

  ost::mol::EntityHandle return_handle = ost::mol::CreateEntity();
  ost::mol::XCSEditor ed = return_handle.EditXCS();
  ost::mol::ChainHandle ch = ed.InsertChain("M");
  ost::mol::ResidueHandle res;

  Real extent = pow(1.0/density,1.0/3.0);

  geom::Vec3 extents = cuboid.GetMax()-cuboid.GetMin();

  int x_atoms = int(extents[0]/extent);
  int y_atoms = int(extents[1]/extent);
  int z_atoms = int(extents[2]/extent);

  geom::Vec3 pos;
  geom::Vec3 origin = grid.GetOrigin();

  geom::Vec3 normalized_a = geom::Vec3(1,0,0);
  geom::Vec3 normalized_b = geom::Vec3(0,1,0);
  geom::Vec3 normalized_c = geom::Vec3(0,0,1);

  geom::Vec3 random_vec;

  timeval t1;
  gettimeofday(&t1, NULL);
  std::srand(t1.tv_usec * t1.tv_usec);


  std::vector<geom::Vec3> positions;

  for(int i=0;i<=x_atoms;++i){
    for(int j=0;j<=y_atoms;++j){
      for(int k=0;k<=z_atoms;++k){
        random_vec[0] = Real(std::rand())/RAND_MAX*extent-extent/2.0;
        random_vec[1] = Real(std::rand())/RAND_MAX*extent-extent/2.0;
        random_vec[2] = Real(std::rand())/RAND_MAX*extent-extent/2.0;
        pos = origin+i*extent*normalized_a+j*extent*normalized_b+k*extent*normalized_c + random_vec;
        if (grid.IsFilled(pos)){
          positions.push_back(pos);
        }
      }
    }
  }


  int actual_number = 0;
  char actual_char='A';

  for(std::vector<geom::Vec3>::iterator it = positions.begin();it!=positions.end();++it){
    res = ed.AppendResidue(ch,"MEM",ost::mol::ResNum(actual_number,actual_char));
    ed.InsertAtom(res,"MEM",*it,"C");

    if(actual_char=='Z'){
      actual_char='A';
      ++actual_number;
    }
    else{
      ++actual_char;
    }

  }

  
  return return_handle;
}

