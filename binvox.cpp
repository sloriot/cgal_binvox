#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/OBJ_reader.h>
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
  Mesh mesh;

  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > triangles;

  std::ifstream in(argv[1]);
  if (!in)
  {
    std::cerr << "ERROR: cannot open argv[1]\n";
    return 1;
  }

  CGAL::read_OBJ(in, points, triangles);

  if ( !PMP::is_polygon_soup_a_polygon_mesh(triangles) )
  {
    std::cerr << "ERROR: input file is not a valid polygon mesh but a soup\n";
    return 1;
  }

  PMP::polygon_soup_to_polygon_mesh(points, triangles, mesh);
  PMP::triangulate_faces(mesh);

  // TODO: extract from command line
  const std::size_t grid_size = 10;
  double gs_d = grid_size;
  CGAL::Bbox_3 bbox = PMP::bbox(mesh);
  // identify the longest side
  double extend = bbox.xmax() - bbox.xmin();
  extend = (std::max)(extend, bbox.ymax() - bbox.ymin());
  extend = (std::max)(extend, bbox.zmax() - bbox.zmin());
  Point_3 origin(bbox.xmin(), bbox.ymin(), bbox.zmin());

  std::ofstream out("output.bin", std::ios_base::binary); // TODO: extract name from commandline
  out << "#binvox 1\n";
  out << "dim " << grid_size << " " << grid_size << " " << grid_size << "\n";
  out << "translate " << origin.x() << " " << origin.y() << " " << origin.z() << "\n";
  out << "scale " << extend << "\n";
  out << "data\n";

  CGAL::Side_of_triangle_mesh<Mesh, K> side_of(mesh);
  
  extend/=gs_d; // TODO rename
  typedef unsigned char byte;
  bool prev=false;
  byte count=0;
  for (double xi=0; xi<gs_d; ++xi)
    for(double zi=0; zi<gs_d; ++zi)
      for(double yi=0; yi<gs_d; ++yi)
      {
        Point_3 query = origin + K::Vector_3((xi+0.5)*extend, (yi+0.5)*extend, (zi+0.5)*extend);
        bool inside = query.x()<=bbox.xmax() && query.y()<=bbox.ymax() && query.z()<=bbox.zmax()
                      && side_of(query) != CGAL::ON_UNBOUNDED_SIDE;
        if (xi==0 && yi==0 && zi==0)
          prev = inside;
        if (prev==inside)
          ++count;
        else
        {
          byte value=prev;
          out << value << count;
          prev=inside;
        }
      }
  byte value=prev;
  out << value << count;

  return 0;
}
