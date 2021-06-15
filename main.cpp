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

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

typedef std::pair<Point, Vector> PointVectorPair;
typedef std::vector<PointVectorPair> PointList;

typedef CGAL::cpp11::array<double, 6> Covariance;

int main()
{
  std::string filename = "/home/jaired/Datasets/PointClouds/ABQ-215-1m-Meru3.ply";
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
  // Estimates normal direction
  const int nb_neighbors = 18;
  CGAL::pca_estimate_normals<CGAL::Sequential_tag>(points.begin(), points.end(),
                  CGAL::First_of_pair_property_map<PointVectorPair>(),
                  CGAL::Second_of_pair_property_map<PointVectorPair>(),
                  nb_neighbors);
  std::cout << "done." << std::endl;
  std::cout << "Orienting normals...";
  // Orients normals
  std::list<PointVectorPair>::iterator unoriented_points_begin =
      CGAL::mst_orient_normals(points.begin(), points.end(),
                  CGAL::First_of_pair_property_map<PointVectorPair>(),
                  CGAL::Second_of_pair_property_map<PointVectorPair>(),
                  nb_neighbors);
  std::cout << "done." << std::endl;
  // Delete unoriented if needed
//  points.erase(unoriented_points_begin, points.end());
  return 0;
}
