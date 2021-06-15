#include <iostream>
#include <utility>
#include <vector>
#include <fstream>

#include <CGAL/IO/read_ply_points.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/vcm_estimate_edges.h>
#include <CGAL/vcm_estimate_normals.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_3.h>
#include <CGAL/compute_average_spacing.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef std::pair<Point, Vector> PointVectorPair;
typedef std::vector<PointVectorPair> PointList;

typedef CGAL::cpp11::array<double, 6> Covariance;

typedef CGAL::Parallel_if_available_tag Concurrency_tag;

int main()
{
  std::string filename = "S:\\Documents\\Graduate School\\REU\\ABQ-215-1m-Meru3.ply";
  std::list<PointVectorPair> points;
  std::ifstream stream(filename);

  std::cout << "Reading file " << filename << "...";
  if (!stream || !CGAL::read_ply_points(stream, std::back_inserter(points),
                          CGAL::First_of_pair_property_map<PointVectorPair>()))
  {
    std::cerr << "Can't read input file." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "done." << std::endl;
  std::cout << "Estimating normal direction...";

  const int nb_neighbors = 18;

  double spacing
      = CGAL::compute_average_spacing<Concurrency_tag>
      (points, nb_neighbors,
          CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()));
  // Then, estimate normals with a fixed radius
  CGAL::pca_estimate_normals<Concurrency_tag>
      (points,
          0, // when using a neighborhood radius, K=0 means no limit on the number of neighbors returns
          CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
          normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()).
          neighbor_radius(2. * spacing)); // use 2*spacing as neighborhood radius
  std::cout << "done." << std::endl;
  std::cout << "Orienting normals...";
 // Orients normals.
    // Note: mst_orient_normals() requires a range of points
    // as well as property maps to access each point's position and normal.
  std::list<PointVectorPair>::iterator unoriented_points_begin =
      CGAL::mst_orient_normals(points, nb_neighbors,
          CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
          normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
  std::cout << "done." << std::endl;
  // Optional: delete points with an unoriented normal
  // if you plan to call a reconstruction algorithm that expects oriented normals.
  //points.erase(unoriented_points_begin, points.end());
  return EXIT_SUCCESS;
}
