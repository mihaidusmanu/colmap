#ifndef COLMAP_SRC_UTIL_FRUSTUM_H_
#define COLMAP_SRC_UTIL_FRUSTUM_H_

#include <vector>

#include <Eigen/Core>

#include "base/camera.h"
#include "base/image.h"

namespace colmap {

typedef Eigen::Matrix<double, 2, 3, Eigen::ColMajor> Matrix2x3d;

class Quadrilateral {
  public:
    Quadrilateral(const Eigen::Vector3d& top_left, const Eigen::Vector3d& top_right, const Eigen::Vector3d& bottom_left, const Eigen::Vector3d& bottom_right);

    const Eigen::Vector3d& TopLeft() const;
    const Eigen::Vector3d& TopRight() const;
    const Eigen::Vector3d& BottomLeft() const;
    const Eigen::Vector3d& BottomRight() const;

    const Eigen::Vector3d Normal() const;

    double Area() const;

    double DistanceToPoint(const Eigen::Vector3d& point) const;

    bool ContainsPointInHalfspace(const Eigen::Vector3d& point) const;

  private:
    const Eigen::Vector3d top_left_;
    const Eigen::Vector3d top_right_;
    const Eigen::Vector3d bottom_left_;
    const Eigen::Vector3d bottom_right_;

    Eigen::Vector3d plane_offset_;
    Matrix2x3d plane_projection_matrix_;
    
    void Check() const;
};

class Pyramid{
  public:
    Pyramid(const Eigen::Vector3d& apex, const Quadrilateral& base);

    double Volume() const;

  private:
    const Eigen::Vector3d apex_;
    const Quadrilateral base_;
};

class Frustum {
  public:
    Frustum(const Eigen::Vector3d& apex, const Quadrilateral& close_face, const Quadrilateral& far_face);

    double Volume();

    bool ContainsPoint(const Eigen::Vector3d& point) const;

    const std::pair<double, double> GetBounds(const size_t coord_idx) const;

  private:
    const Eigen::Vector3d apex_;
    const Quadrilateral close_face_;
    const Quadrilateral far_face_;
    
    double volume_;
};

}  // namespace colmap

#endif  // COLMAP_SRC_UTIL_FRUSTUM_H_
