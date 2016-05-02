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

/*
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
*/

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

    min_x_idx = int(std::floor((min_x-origin_[0])/bin_size_));
    max_x_idx = int(std::ceil((max_x-origin_[0])/bin_size_));
    min_y_idx = int(std::floor((min_y-origin_[1])/bin_size_));
    max_y_idx = int(std::ceil((max_y-origin_[1])/bin_size_));
    min_z_idx = int(std::floor((min_z-origin_[2])/bin_size_));
    max_z_idx = int(std::ceil((max_z-origin_[2])/bin_size_));

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
    throw ost::Error("Cannot calculate solvated surface in grid without attached surface!");
  }

  if(!is_flooded_){
    this->Flood();
    is_flooded_=true;
  }

  solvated_surf_ = ost::mol::CreateSurface();
  solvated_vertices_.clear();

  int level,x,y;

  /*
  for(level=0;level<zbins_;++level){
    for(x=2;x<xbins_-2;++x){
      for(y=2;y<ybins_-2;++y){
        if(grid_[x][y][level] == 1){
          //we check the first layer
          if(grid_[x-1][y][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x+1][y][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x][y-1][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x][y+1][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x+1][y+1][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x+1][y-1][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x-1][y+1][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x-1][y-1][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          //now we check the second layer
          if(grid_[x-2][y-2][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x-2][y-1][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }          
          if(grid_[x-2][y][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x-2][y+1][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }          
          if(grid_[x-2][y+2][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x-1][y-2][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }          
          if(grid_[x-1][y+2][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x][y-2][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }          
          if(grid_[x][y+2][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }          
          if(grid_[x+1][y-2][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x+1][y+2][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }          
          if(grid_[x+2][y-2][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x+2][y-1][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }          
          if(grid_[x+2][y][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x+2][y+1][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
          if(grid_[x+2][y+2][level] == 2){
            grid_[x][y][level] = 10;
            continue;
          }
        }
      }
    }
  }
  */


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

    min_x_idx = int(std::floor((min_x-origin_[0])/bin_size_));
    max_x_idx = int(std::ceil((max_x-origin_[0])/bin_size_));
    min_y_idx = int(std::floor((min_y-origin_[1])/bin_size_));
    max_y_idx = int(std::ceil((max_y-origin_[1])/bin_size_));
    min_z_idx = int(std::floor((min_z-origin_[2])/bin_size_));
    max_z_idx = int(std::ceil((max_z-origin_[2])/bin_size_));

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
              solvated_vertices_.push_back(tri.v0);
              solvated_vertices_.push_back(tri.v1);
              solvated_vertices_.push_back(tri.v2);
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


ost::mol::SurfaceHandle TransformSurface(ost::mol::SurfaceHandle& s, geom::Transform& t){
  
  
  ost::mol::SurfaceHandle transformed_surf = ost::mol::CreateSurface();
  ost::mol::SurfaceTriIDList tri_id_list = s.GetTriIDList();
  ost::mol::SurfaceVertexIDList v_id_list = s.GetVertexIDList();

  std::map<ost::mol::SurfaceVertexID, ost::mol::SurfaceVertexID> v_mapper;

  ost::mol::SurfaceVertex v;
  ost::mol::SurfaceVertex transformed_v;
  geom::Vec3 transformed_pos;
  geom::Vec3 transformed_norm;
  geom::Vec3 norm_helper;
  ost::mol::SurfaceTri tri;


  //vertices get added first
  for(ost::mol::SurfaceVertexIDList::iterator i = v_id_list.begin(); i!=v_id_list.end(); ++i){
    v = s.GetVertex(*i);
    norm_helper = t.Apply(v.position+v.normal);
    transformed_pos = t.Apply(v.position);
    transformed_norm = norm_helper-transformed_pos;
    transformed_v = ost::mol::SurfaceVertex(transformed_pos,transformed_norm,
                                            v.type, v.atom);
    v_mapper[*i] = transformed_surf.AddVertex(transformed_v);
  }

  //let's add the faces...
  for(ost::mol::SurfaceTriIDList::iterator i=tri_id_list.begin();i!=tri_id_list.end();++i){ 
    tri = s.GetTri(*i);
    transformed_surf.AddTri(v_mapper[tri.v0],v_mapper[tri.v1],v_mapper[tri.v2]);
  }
  
  return transformed_surf;
}


geom::Vec3 RotateAroundAxis(geom::Vec3 point, geom::Vec3 axis, Real angle){

  Real aa,ab,ac,ba,bb,bc,ca,cb,cc,one_m_cos,cos_ang,sin_ang;
  geom::Mat3 rot;
  geom::Vec3 result;

  cos_ang = std::cos(angle);
  sin_ang = std::sin(angle);
  one_m_cos = 1-cos_ang;

  aa = cos_ang+axis[0]*axis[0]*one_m_cos;
  ab = axis[0]*axis[1]*one_m_cos-axis[2]*sin_ang;
  ac = axis[0]*axis[2]*one_m_cos+axis[1]*sin_ang;

  ba = axis[1]*axis[0]*one_m_cos+axis[2]*sin_ang;
  bb = cos_ang+axis[1]*axis[1]*one_m_cos;
  bc = axis[1]*axis[2]*one_m_cos-axis[0]*sin_ang;

  ca = axis[2]*axis[0]*one_m_cos-axis[1]*sin_ang;
  cb = axis[2]*axis[1]*one_m_cos+axis[0]*sin_ang;
  cc = cos_ang+axis[2]*axis[2]*one_m_cos;

  rot = geom::Mat3(aa,ab,ac,ba,bb,bc,ca,cb,cc);

  result = rot*point;

  return result;

}

geom::Vec3 GetTiltAxis(geom::Vec3& x){
  
  geom::Vec3 base = geom::Normalize(x);
  geom::Vec3 tilt_axis;

  if(std::abs(base[0])>0.9999){
    tilt_axis = geom::Vec3(0.0,1.0,0.0);
  }
  else{
    geom::Vec3 x_base = geom::Vec3(1.0,0.0,0.0);
    tilt_axis = geom::Normalize(geom::Cross(base,x_base));
  }

  return tilt_axis;
}

FindMemParam MinimizeAlongZ(std::vector<geom::Vec3>& atom_positions, std::vector<Real>& transfer_energies, Real lambda){


  geom::Vec3 normalized_axis(0.0,0.0,1.0);
  geom::Vec3 tilt_axis(1.0,0.0,0.0);

  //top ten of the grid search solutions will be saved in here
  std::vector<FindMemParam> initial_solutions;
  std::pair<std::pair<Real,Real>, Real> actual_solution;

  //the minimal solution will be stored in here
  FindMemParam final_solution;
  final_solution.energy = std::numeric_limits<Real>::max();

  geom::Vec3 tilted_axis;
  geom::Vec3 scan_axis;
  Real tilt_rad,angle_rad;

  FindMemParam parameters;

  for(int tilt_deg = 0; tilt_deg<=52; tilt_deg+=4){

    tilt_rad = Real(tilt_deg)/360*2*M_PI;
    tilted_axis = RotateAroundAxis(normalized_axis,tilt_axis,tilt_rad);

    //if we would choos equidistant angles, rotations from higher tilts would be sampled more sparse.
    //We tilt until 52 degrees. If we fully rotate with a tilted axis in a sphere of radius 1, we would
    //get a circle with circumference of pi. Let's say we want to sample this region widh d_angle of
    //4 degrees, we need 90 different angles and the points are on average 
    //0.0550135287662 (arclength) apart. We can now calculate the circumference for every circle at
    //a specific tilt and then have a look how many angles we need to get a similar d_arclength...

    Real circumference = 2*M_PI*std::sin(tilt_rad);
    int num_angles = int(ceil(std::max(12,int(ceil(circumference/0.0550135287662)))));
    Real d_angle = 2*M_PI/num_angles;

    for(int actual_angle = 0; actual_angle<num_angles; ++actual_angle){
      
      angle_rad = actual_angle*d_angle;
      scan_axis = RotateAroundAxis(tilted_axis,normalized_axis,angle_rad);

      actual_solution = ScanAxis(atom_positions,transfer_energies,scan_axis);

      if(initial_solutions.size()==0){
        parameters.axis = normalized_axis;
        parameters.tilt_axis = tilt_axis;
        parameters.tilt = tilt_rad;
        parameters.angle = angle_rad;
        parameters.width = actual_solution.first.first;
        parameters.pos = actual_solution.first.second;
        parameters.energy = actual_solution.second;
        initial_solutions.push_back(parameters);
        //happens only once... At this time the tilt is zero
        //and it makes no sense to do additional angle rotations...
        break;
      }

      if(initial_solutions.size()<10){

        parameters.axis = normalized_axis;
        parameters.tilt_axis = tilt_axis;
        parameters.tilt = tilt_rad;
        parameters.angle = angle_rad;
        parameters.width = actual_solution.first.first;
        parameters.pos = actual_solution.first.second;
        parameters.energy = actual_solution.second;

        //we either append or insert the parameters
        if(actual_solution.second>=initial_solutions.back().energy){
          initial_solutions.push_back(parameters);
          continue;
        }

        for(std::vector<FindMemParam>::iterator i = initial_solutions.begin();
            i!=initial_solutions.end();++i){
          if(i->energy > actual_solution.second){
            initial_solutions.insert(i,parameters);
            break;
          }
        }
        continue;
      }

      if(actual_solution.second < initial_solutions.back().energy){
        parameters.axis = normalized_axis;
        parameters.tilt_axis = tilt_axis;
        parameters.tilt = tilt_rad;
        parameters.angle = angle_rad;
        parameters.width = actual_solution.first.first;
        parameters.pos = actual_solution.first.second;
        parameters.energy = actual_solution.second;
        for(std::vector<FindMemParam>::iterator i= initial_solutions.begin();
            i!=initial_solutions.end();++i){
          if(i->energy > actual_solution.second){
            initial_solutions.insert(i,parameters);
            initial_solutions.pop_back();
            break;
          }
        }
      } 
    }
  }

  EnergyF en_f(atom_positions, transfer_energies, 0.9, 0.0);

  Eigen::Matrix<Real, 4, 1> lm_parameters;
  FindMemParam temp;
  ost::img::alg::LevenbergMarquardt<EnergyF,EnergyDF>::Results lm_result; 

  //go through the initial solutions found by simple scanning and run
  //levenberg marquardt minimization
  for(int i=0;i<5;++i){

    temp = initial_solutions[i];

    lm_parameters(0,0) = temp.tilt;
    lm_parameters(1,0) = temp.angle;
    lm_parameters(2,0) = temp.width;
    lm_parameters(3,0) = temp.pos;

    //set offset parameter
    en_f.offset = std::abs(temp.energy*2);

    ost::img::alg::LevenbergMarquardt<EnergyF,EnergyDF> lm(en_f);
    lm_result = lm.minimize(&lm_parameters);
    Real minimized_energy = en_f(lm_parameters)(0,0)-en_f.offset;
    if(minimized_energy<final_solution.energy){
      FindMemParam mem_param;
      mem_param.tilt = lm_parameters(0,0);
      mem_param.angle = lm_parameters(1,0);
      mem_param.width = lm_parameters(2,0);
      mem_param.pos = lm_parameters(3,0);
      mem_param.energy = minimized_energy;
      mem_param.axis = temp.axis;
      mem_param.tilt_axis = temp.tilt_axis;
      final_solution = mem_param;
    }
  }
  return final_solution;
}

std::pair<std::pair<Real,Real>, Real> ScanAxis(std::vector<geom::Vec3>& atom_positions, std::vector<Real>& transfer_energies, geom::Vec3& axis){
  
  geom::Vec3 normalized_axis = geom::Normalize(axis);
  std::vector<Real> pos_on_axis;

  for(std::vector<geom::Vec3>::iterator i = atom_positions.begin();
      i != atom_positions.end(); ++i){
    pos_on_axis.push_back(geom::Dot(*i,normalized_axis));
  }

  Real min_pos = pos_on_axis[0];
  Real max_pos = min_pos;

  for(std::vector<Real>::iterator i = pos_on_axis.begin(); 
      i != pos_on_axis.end(); ++i){
    if(min_pos>(*i)) min_pos = (*i);
    if(max_pos<(*i)) max_pos = (*i);
  }

  int width = int(ceil(max_pos)-floor(min_pos));


  //energies representing the energy profile along the axis
  Real energies[width];

  for(int i=0;i<width;++i){
    energies[i]=0.0;
  }

  std::vector<Real>::iterator i,j;

  for(i = pos_on_axis.begin(),j = transfer_energies.begin();
      i != pos_on_axis.end() && j != transfer_energies.end(); ++i, ++j){
    energies[int(floor((*i)-min_pos))] += (*j);
  }

  // int pair: position and width, real value: energy
  std::pair<std::pair<Real,Real>, Real> solution = std::make_pair(std::make_pair(0.0,0.0),std::numeric_limits<Real>::max());
  Real energy;

  for(int window_width = 10; window_width<=40; window_width+=1){

    if(window_width>width) break;

    energy=0.0;
    for(int i=0;i<window_width;++i){
      energy+=energies[i];
    }
    if(energy<solution.second){
      Real center = min_pos+Real(window_width)/2;
      solution = std::make_pair(std::make_pair(Real(window_width),center),energy);
    }
    
    for(int pos = 1; pos<=width-window_width; ++pos){
      energy-=energies[pos-1];
      energy+=energies[pos+window_width-1];
      if(energy<solution.second){
        Real center = min_pos+pos+Real(window_width)/2;
        solution = std::make_pair(std::make_pair(Real(window_width),center),energy);
      }
    }
  }
  return solution;
}


FindMemParam FindMembrane(ost::mol::EntityHandle& ent, ost::mol::SurfaceHandle& surf, std::vector<Real>& asa){

  if(ent.GetAtomCount() != static_cast<int>(asa.size())){
    throw ost::io::IOException("Number of atoms and number of elements in asa must be consistent!");
  }

  ost::io::IOProfile io_profile;


  ost::mol::EntityHandle rotated_ent;
  ost::mol::SurfaceHandle rotated_surf;

  std::vector<Real> transfer_energies;
  std::vector<geom::Vec3> atom_positions;

  FindMemParam final_solution;
  FindMemParam actual_solution;

  final_solution.energy = std::numeric_limits<Real>::max();

  //The whole grid stuff for defining surface atoms only works along the
  //z-axis. We therefore have to transform the structure. We use a rotation 
  //around the z-axis with subsequent rotation around the x-axis for this task
  std::vector<Real> z_rotations;
  std::vector<Real> x_rotations;

  z_rotations.push_back(0.0);
  z_rotations.push_back(0.0);
  z_rotations.push_back(1*(2*M_PI/5));
  z_rotations.push_back(2*(2*M_PI/5));
  z_rotations.push_back(3*(2*M_PI/5));
  z_rotations.push_back(4*(2*M_PI/5));
  
  x_rotations.push_back(0.0);
  x_rotations.push_back(M_PI/4);
  x_rotations.push_back(M_PI/4);
  x_rotations.push_back(M_PI/4);
  x_rotations.push_back(M_PI/4);
  x_rotations.push_back(M_PI/4);

  
  geom::Transform transform;
  geom::Mat3 z_rot_matrix, x_rot_matrix, rot_matrix;
  Real z_rot, x_rot;


  //we look, which vertex is associated to which atom...
  //this will speed up some calculations later on
  std::map<ost::mol::SurfaceVertexID,int> vertex_atom_mapper;

  ost::mol::SurfaceVertexIDList v_list = surf.GetVertexIDList();
  ost::mol::SurfaceVertex vert;
  ost::mol::AtomHandleList a_list;

  for(ost::mol::SurfaceVertexIDList::iterator it=v_list.begin();it!=v_list.end();++it){
    vert = surf.GetVertex(*it);
    a_list = ent.FindWithin(vert.position,3.0);
    if(a_list.empty()) continue;
    ost::mol::AtomHandleList::iterator ait=a_list.begin();
    Real best_dist2 = geom::Length2(ait->GetPos()-vert.position);
    ost::mol::AtomHandle best_atom = *ait;
    ++ait;
    for(;ait!=a_list.end();++ait) {
      Real dist2 = geom::Length2(ait->GetPos()-vert.position);
      if(dist2<best_dist2) {
        best_dist2=dist2;
        best_atom = *ait;
      }
    }
    best_atom.SetIntProp("surface_vertex",*it);
  }

  a_list = ent.GetAtomList();

  for(int i = 0; i < a_list.size(); ++i){
    if(a_list[i].HasProp("surface_vertex")){
      vertex_atom_mapper[a_list[i].GetIntProp("surface_vertex")] = i;
    }
  }


  for(int t_counter = 0; t_counter<z_rotations.size();++t_counter){

    //rotate entity and surface to start a search around the global z-axis
    rotated_ent = ent.Copy();
    ost::mol::XCSEditor ed = rotated_ent.EditXCS();

    z_rot = z_rotations[t_counter];
    x_rot = x_rotations[t_counter];

    z_rot_matrix = geom::Mat3(std::cos(z_rot),-std::sin(z_rot),0.0,
                              std::sin(z_rot),std::cos(z_rot),0.0,
                              0.0, 0.0, 1.0);

    x_rot_matrix = geom::Mat3(1.0, 0.0, 0.0,
                              0.0, std::cos(x_rot),-std::sin(x_rot),
                              0.0,std::sin(x_rot),std::cos(x_rot));

    //rot_matrix = z_rot_matrix*x_rot_matrix;
    rot_matrix = x_rot_matrix*z_rot_matrix;
    transform.SetRot(rot_matrix);

    ed.ApplyTransform(transform.GetMatrix());
    rotated_surf = TransformSurface(surf, transform);

    //to define the grid size later on, we have to find the structures extent
    ost::mol::AtomHandleList a_list = rotated_ent.GetAtomList();
    Real max_x = std::numeric_limits<Real>::min();
    Real min_x = std::numeric_limits<Real>::max();
    Real max_y = std::numeric_limits<Real>::min();
    Real min_y = std::numeric_limits<Real>::max();
    Real max_z = std::numeric_limits<Real>::min();
    Real min_z = std::numeric_limits<Real>::max();
    geom::Vec3 pos;
    for(ost::mol::AtomHandleList::iterator i = a_list.begin(); i!=a_list.end();++i){
      pos = i->GetPos();
      if(pos[0]>max_x) max_x = pos[0];
      if(pos[0]<min_x) min_x = pos[0];
      if(pos[1]>max_y) max_y = pos[1];
      if(pos[1]<min_y) min_y = pos[1];
      if(pos[2]>max_z) max_z = pos[2];
      if(pos[2]<min_z) min_z = pos[2];
    }

    max_x += 2;
    min_x -= 2;
    max_y += 2;
    min_y -= 2;
    max_z += 2;
    min_z -= 2;

    //define a cuboid with the previously found extents and generate a grid
    //with that size
    geom::Vec3 min_vec(min_x,min_y,min_z);
    geom::Vec3 max_vec(max_x,max_y,max_z);
    geom::AlignedCuboid cuboid(min_vec,max_vec);
    SolvationGrid grid(cuboid, 0.5);

    //Add the rotated surface to the grid
    grid.AddSurface(rotated_surf);
    grid.Flood();
    grid.CalculateSolvatedSurface();

    //note, that we defined the mapping between the vertex indices and the
    //corresponding atoms. It is therefore possible, to just look at the
    //remaining vertices of the solvated surface
    v_list = grid.GetSolvatedVertices();
    a_list = rotated_ent.GetAtomList();

    for(ost::mol::SurfaceVertexIDList::iterator i = v_list.begin();
        i != v_list.end(); ++i){
      a_list[vertex_atom_mapper[*i]].SetBoolProp("surface",true);
    }

    /*
    std::stringstream ss;
    ss  << "test" << t_counter << ".pdb";

    ost::mol::EntityView solvated_atoms = rotated_ent.Select("gasurface:false=true");

    ost::io::PDBWriter writer(ss.str(), io_profile);

    writer.Write(solvated_atoms);
    */


    String element;
    unsigned char bond_order;
    ost::mol::BondHandleList bond_list;
    bool assigned_energy=false;

    ost::mol::AtomHandleList::iterator i;
    std::vector<Real>::iterator j;

    transfer_energies.clear();
    atom_positions.clear();

    for(i=a_list.begin(), j=asa.begin();i!=a_list.end(),j!=asa.end();++i,++j){
      //we don't even have to check for true, since prop is only set for surface atoms
      if(!i->HasProp("surface")){
        continue;
      }
      atom_positions.push_back(i->GetPos());
      element = i->GetElement();
      if(element=="S") transfer_energies.push_back((*j)*(10.0));
      else if(element=="N") transfer_energies.push_back((*j)*(53.0));
      else if(element=="O") transfer_energies.push_back((*j)*(57.0));
      else if(element=="C"){
        assigned_energy=false;
        bond_list = i->GetBondList();
        for(ost::mol::BondHandleList::iterator k=bond_list.begin();k!=bond_list.end();++k){
          bond_order = k->GetBondOrder();
          if(bond_order>'1'){
            transfer_energies.push_back((*j)*(-19.0));
            assigned_energy=true;
            break;;
          }
        }
        if(!assigned_energy) transfer_energies.push_back((*j)*(-22.6));
      }
      else if(element=="H") continue;
      else transfer_energies.push_back(0.0); // in case of unknown atom...
    }

    //we now try to find a solution, taking the global z-axis as a starting point
    actual_solution = MinimizeAlongZ(atom_positions,transfer_energies,0.9);

    //Check whether newly calculated minimum is better than global solution
    if(actual_solution.energy < final_solution.energy){
      final_solution = actual_solution;
      //we still have to correct the membrane and tilt axis for the initial transform
      final_solution.tilt_axis = transform.ApplyInverse(final_solution.tilt_axis);
      final_solution.axis = transform.ApplyInverse(final_solution.axis);

    }
  }
  return final_solution;
}

EnergyF::FMatrixType EnergyF::operator()(const Eigen::Matrix<Real, 4, 1>& x) const{

  FMatrixType result;
  result(0,0) = 0.0;
  geom::Vec3 tilted_axis = axis;

  tilted_axis = RotateAroundAxis(tilted_axis, tilt_axis, x(0,0));
  tilted_axis = RotateAroundAxis(tilted_axis, axis, x(1,0));
  

  std::vector<Real>::const_iterator e_it = transfer_energies.begin();
  std::vector<geom::Vec3>::const_iterator pos_it = positions.begin();

  Real distance_to_center;
  Real pos_on_axis;
  Real half_width = x(2,0)/2.0;
  Real exponent;
  Real energy;

  for(;pos_it!=positions.end() && e_it!=transfer_energies.end();++pos_it,++e_it){
    pos_on_axis = geom::Dot(tilted_axis,*pos_it);
    distance_to_center = std::abs(x[3]-pos_on_axis);
    exponent = (distance_to_center-half_width)/lambda;
    energy = (1.0/(1.0+std::exp(exponent)))*(*e_it);
    result(0,0)=result(0,0)+energy;
  }

  result(0,0)+=offset;

  return result;
}

Eigen::Matrix<Real,1,4> EnergyDF::operator()(const Eigen::Matrix<Real, 4, 1>& x) const{

  Eigen::Matrix<Real,1,4> result;
  Eigen::Matrix<Real, 4, 1> parameter1 = x;
  Eigen::Matrix<Real, 4, 1> parameter2 = x;

  parameter1(0,0)+=d_tilt;
  parameter2(0,0)-=d_tilt;
  result(0,0) = (function(parameter1)(0,0)-function(parameter2)(0,0))/(2*d_tilt);

  parameter1=x;
  parameter2=x;
  parameter1(1,0)+=d_angle;
  parameter2(1,0)-=d_angle;
  result(0,1) = (function(parameter1)(0,0)-function(parameter2)(0,0))/(2*d_angle);

  parameter1=x;
  parameter2=x;
  parameter1(2,0)+=d_width;
  parameter2(2,0)-=d_width;
  result(0,2) = (function(parameter1)(0,0)-function(parameter2)(0,0))/(2*d_width);

  parameter1=x;
  parameter2=x;
  parameter1(3,0)+=d_pos;
  parameter2(3,0)-=d_pos;
  result(0,3) = (function(parameter1)(0,0)-function(parameter2)(0,0))/(2*d_pos);
  
  return result;
}


geom::Vec3 FindMemParam::GetMembraneAxis(){  

  geom::Vec3 result = RotateAroundAxis(axis,tilt_axis,tilt);
  result = RotateAroundAxis(result,axis,angle);
  return result;

}

FindMemParam FindMemParam::Load(const String& filename){

  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open MemParam. File '"
       << filename << "' does not exist";
    throw ost::io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  ost::io::BinaryDataSource ds(stream);
  FindMemParam param;
  ds >> param;
  return param;
}


void FindMemParam::Save(const String& filename){
  std::ofstream stream(filename.c_str(), std::ios_base::binary);
  std::stringstream ss;
  if(! stream){
    ss<<"Could not write to file '";
    ss<<filename<<"'";
    throw ost::io::IOException(ss.str());
  } 
  ost::io::BinaryDataSink ds(stream);
  ds << *this;
}




