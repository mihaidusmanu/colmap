#include "util/frustum.h"

// TODO: All this needs to be tested :)!

namespace colmap {

double epsilon = 1e-6;

Quadrilateral::Quadrilateral(
  const Eigen::Vector3d& top_left,
  const Eigen::Vector3d& top_right,
  const Eigen::Vector3d& bottom_left,
  const Eigen::Vector3d& bottom_right
) : top_left_(top_left), top_right_(top_right), bottom_left_(bottom_left), bottom_right_(bottom_right) {
  // TODO: Corner names might be confusing... Change them?

  Check();

  // Plane offset.
  plane_offset_ = top_left_;

  // Compute orthonormal plane projection matrix.
  // First basis vector.
  Eigen::Vector3d v1 = top_right_ - top_left_;
  v1.normalize();

  // Second basis vector.
  Eigen::Vector3d v2 = bottom_left_ - top_left_;
  v2 -= v2.dot(v1) * v1;
  v2.normalize();

  Matrix2x3d P;
  P << v1.transpose(), v2.transpose();

  plane_projection_matrix_ = P;
}

void Quadrilateral::Check() const {
  Eigen::Matrix3d M;
  M << top_right_ - top_left_, bottom_left_ - top_left_, bottom_right_ - top_left_;
  CHECK_LE(M.determinant(), epsilon);
}

const Eigen::Vector3d& Quadrilateral::TopLeft() const {
  return top_left_;
}

const Eigen::Vector3d& Quadrilateral::TopRight() const {
  return top_right_;
}

const Eigen::Vector3d& Quadrilateral::BottomLeft() const {
  return bottom_left_;
}

const Eigen::Vector3d& Quadrilateral::BottomRight() const {
  return bottom_right_;
}

const Eigen::Vector3d Quadrilateral::Normal() const {
  return plane_projection_matrix_.row(0).cross(plane_projection_matrix_.row(1));
}

double Quadrilateral::Area() const {
  // Project to the plane (2D).
  const Eigen::Vector2d top_left = plane_projection_matrix_ * (top_left_ - plane_offset_);
  const Eigen::Vector2d top_right = plane_projection_matrix_ * (top_right_ - plane_offset_);
  const Eigen::Vector2d bottom_left = plane_projection_matrix_ * (bottom_left_ - plane_offset_);
  const Eigen::Vector2d bottom_right = plane_projection_matrix_ * (bottom_right_ - plane_offset_);

  // A = 1 / 2 * p * q * sin(theta)
  // p, q  - lengths of diagonals
  // theta - angle between diagonals
  
  // Compute theta.
  const Eigen::Vector3d diag1 = top_left.homogeneous().cross(bottom_right.homogeneous());
  const Eigen::Vector3d diag2 = top_right.homogeneous().cross(bottom_left.homogeneous());
  // Intersection of diagonals.
  const Eigen::Vector2d inters = diag1.cross(diag2).hnormalized();

  // theta = <(top_left, inters, top_right)
  // Law of cosines: c^2 = a^2 + b^2 - 2 a b cos(theta)
  const double a = (top_left - inters).norm();
  const double b = (top_right - inters).norm();
  const double c = (top_left - top_right).norm();
  const double theta = acos((a * a + b * b - c * c) / (2.0 * a * b));

  // Compute lengths of diagonals.
  const double p = (top_left - bottom_right).norm();
  const double q = (top_right - bottom_left).norm();

  return 1.0 / 2.0 * p * q * sin(theta);
}

double Quadrilateral::DistanceToPoint(const Eigen::Vector3d& point) const {
  const Eigen::Vector3d point_t = point - plane_offset_;
  const Eigen::Vector3d proj = 
    point_t.dot(plane_projection_matrix_.row(0)) * plane_projection_matrix_.row(0) +
    point_t.dot(plane_projection_matrix_.row(1)) * plane_projection_matrix_.row(1);

  return (point_t - proj).norm();
}

bool Quadrilateral::ContainsPointInHalfspace(const Eigen::Vector3d& point) const {
  const Eigen::Vector3d point_t = point - plane_offset_;
  return (Normal().dot(point_t) >= -epsilon);
}

Pyramid::Pyramid(const Eigen::Vector3d& apex, const Quadrilateral& base) : apex_(apex), base_(base) { };

double Pyramid::Volume() const {
  // V = 1 / 3 * height * base_area.
  return 1.0 / 3.0 * base_.DistanceToPoint(apex_) * base_.Area();
}


Frustum::Frustum(const Eigen::Vector3d& apex, const Quadrilateral& close_face, const Quadrilateral& far_face) : apex_(apex), close_face_(close_face), far_face_(far_face) {
  volume_ = Pyramid(apex_, far_face_).Volume() - Pyramid(apex_, close_face_).Volume();
}

double Frustum::Volume() const {
  return volume_;
}

bool Frustum::ContainsPoint(const Eigen::Vector3d& point) const {
  // Corners should be arranged such that 
  // n = (top_right_ - top_left_) x (bottom_left_ - top_left_)
  // points towards the inside of the frustum.
  return (
    close_face_.ContainsPointInHalfspace(point) &
    Quadrilateral(far_face_.TopRight(), far_face_.TopLeft(), far_face_.BottomRight(), far_face_.BottomLeft()).ContainsPointInHalfspace(point) &
    Quadrilateral(far_face_.TopLeft(), close_face_.TopLeft(), far_face_.BottomLeft(), close_face_.BottomLeft()).ContainsPointInHalfspace(point) &
    Quadrilateral(far_face_.TopRight(), close_face_.TopRight(), far_face_.TopLeft(), close_face_.TopLeft()).ContainsPointInHalfspace(point) &
    Quadrilateral(close_face_.TopRight(), far_face_.TopRight(), close_face_.BottomRight(), far_face_.BottomRight()).ContainsPointInHalfspace(point) &
    Quadrilateral(close_face_.BottomRight(), far_face_.BottomRight(), close_face_.BottomLeft(), far_face_.BottomLeft()).ContainsPointInHalfspace(point));
}

const std::pair<double, double> Frustum::GetBounds(const size_t coord_idx) const {
  CHECK_LE(coord_idx, 2);

  const std::vector<double> values = {
    close_face_.TopLeft()(coord_idx), far_face_.TopLeft()(coord_idx),
    close_face_.TopRight()(coord_idx), far_face_.TopRight()(coord_idx),
    close_face_.BottomLeft()(coord_idx), far_face_.BottomLeft()(coord_idx),
    close_face_.BottomRight()(coord_idx), far_face_.BottomRight()(coord_idx)};
  
  const double min = *std::min_element(values.begin(), values.end());
  const double max = *std::max_element(values.begin(), values.end());

  return std::make_pair(min, max);
}

}  // namespace colmap
