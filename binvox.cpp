#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Side_of_triangle_mesh.h>


#include <fstream>

//QString command("\"" + binvoxProgramFileInfo.absoluteFilePath()+ "\" -pb -d "+ QString::number(voxelizationResolution) + " \"" +scaledFilePath + "\"");

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  // parse options
  int oi=0;
  int options_found=0;
  std::string fname="default.obj";
  std::size_t grid_size=0;
  while(++oi!=argc)
  {
    std::string opt=argv[oi];
    if (opt[0]=='-')
    {
      if(opt=="-d")
      {
        ++oi;
        if (oi==argc)
        {
          std::cerr << "ERROR: no size to -d\n";
          return 1;
        }

        int gs = std::stoi(argv[oi]);
        if (gs < 0 || gs > 1024)
        {
          std::cerr << "ERROR: invalid size given to -d\n";
          return 1;
        }
        grid_size=gs;
        ++options_found;
      }
    }
    else
    {
      // should be the filename
      fname = opt;
      ++options_found;
    }
  }

  if (options_found!=2)
  {
    std::cerr << "ERROR: missing -d SIZE and/or input obj file not provided\n";
    return 1;
  }

  Mesh mesh;

  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > triangles;
 
  if (! CGAL::IO::read_polygon_soup(fname, points, triangles))
  {
     std::cerr << "ERROR: cannot read " << argv[1]  << std::endl;
      return 1;
  }

  if ( !PMP::is_polygon_soup_a_polygon_mesh(triangles) )
  {
    std::cerr << "ERROR: input file is not a valid polygon mesh but a soup\n";
    return 1;
  }

  PMP::polygon_soup_to_polygon_mesh(points, triangles, mesh);
  PMP::triangulate_faces(mesh);

  double gs_d = grid_size;
  CGAL::Bbox_3 bbox = PMP::bbox(mesh);
  // identify the longest side
  double extend = bbox.xmax() - bbox.xmin();
  extend = (std::max)(extend, bbox.ymax() - bbox.ymin());
  extend = (std::max)(extend, bbox.zmax() - bbox.zmin());
  Point_3 origin(bbox.xmin(), bbox.ymin(), bbox.zmin());

  std::string prefix = fname.substr(0, fname.size()-4);
  std::ofstream out(prefix+std::string(".binvox"), std::ios_base::binary);
  out << "#binvox 1\n";
  out << "dim " << grid_size << " " << grid_size << " " << grid_size << "\n";
  out << "translate " << origin.x() << " " << origin.y() << " " << origin.z() << "\n";
  out << "scale " << extend << "\n";
  out << "data\n";

  CGAL::Side_of_triangle_mesh<Mesh, K> side_of(mesh);

  double normed_extend=extend/gs_d;
  typedef unsigned char byte;
  bool prev=false;
  byte count=0;
  for (double xi=0; xi<gs_d; ++xi)
    for(double zi=0; zi<gs_d; ++zi)
      for(double yi=0; yi<gs_d; ++yi)
      {
        Point_3 query = origin + K::Vector_3((xi+0.5)*normed_extend, (yi+0.5)*normed_extend, (zi+0.5)*normed_extend);
        bool inside = query.x()<=bbox.xmax() && query.y()<=bbox.ymax() && query.z()<=bbox.zmax()
                      && side_of(query) != CGAL::ON_UNBOUNDED_SIDE;
        if (xi==0 && yi==0 && zi==0)
          prev = inside;
        if (prev==inside)
        {
          ++count;
          if (count==0) // if too many indentical values, byte type might go around
          {
            --count;
            byte value=prev;
            out << value << count;
            count=1;
          }
        }
        else
        {
          byte value=prev;
          out << value << count;
          prev=inside;
          count=1;
        }
      }
  byte value=prev;
  out << value << count;

  return 0;
}
