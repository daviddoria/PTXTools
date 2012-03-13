#include <iostream>

namespace CameraCalibration
{
  
template<typename T>
T Centroid(const typename std::vector<T,typename Eigen::aligned_allocator<T> >& points)
{
  T centroid;
  centroid.fill(0.0f);

  for(unsigned int i = 0; i < points.size(); ++i)
    {
    centroid += points[i];
    }

  float numberOfPoints = static_cast<float>(points.size());
  centroid /= numberOfPoints;

  return centroid;
}

template<typename T>
Eigen::MatrixXd ComputeNormalizationTransform(const typename std::vector<T,typename Eigen::aligned_allocator<T> >& points)
{
  unsigned int numberOfPoints = points.size();
  unsigned int dimensions = points[0].rows();
  // std::cout << "dimensions: " << dimensions << std::endl;
  
  T centroid = Centroid<T>(points);

  typename std::vector<T,typename Eigen::aligned_allocator<T> > centeredPoints(numberOfPoints);

  // Shift the origin of the points to the centroid
  for(unsigned int i = 0; i < numberOfPoints; ++i)
    {
    centeredPoints[i] = points[i] - centroid;
    }

  // Normalize the points so that the average distance from the origin is equal to sqrt(2).
  // Compute the average distance
  float totalDistance = 0.0f;
  for(unsigned int i = 0; i < numberOfPoints; ++i)
    {
    totalDistance += centeredPoints[i].norm();
    }

  float averageDistance = totalDistance / static_cast<float>(numberOfPoints);

  float scale = sqrt(static_cast<float>(dimensions))/averageDistance;

  Eigen::MatrixXd similarityTransform(dimensions+1, dimensions+1);
  similarityTransform.fill(0.0f);
  similarityTransform(dimensions, dimensions) = 1; // Set the lower right corner to 1

  // Set the diagonal to 'scale'
  for(unsigned int i = 0; i < dimensions; ++i)
    {
    similarityTransform(i,i) = scale;
    }

  // Set the first 'dimensions' entries of the last column to -scale * centroid(row)
  for(unsigned int i = 0; i < dimensions; ++i)
    {
    similarityTransform(i,dimensions) = -scale * centroid(i);
    }

  //std::cout << "SimilarityTransform: " << similarityTransform << std::endl;

  return similarityTransform;
}

} // end namespace
