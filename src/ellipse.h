#pragma once

#include <Eigen/Dense>

class Ellipse_t
{
public:
  /// the center of the ellipse
  Eigen::Vector2d center = {0, 0};
  /// the axes of the ellipse
  Eigen::Vector2d axes = {0, 0};
  /// the rotation of the ellipse from x axis to principal axis
  double rotation = 0;
  // residue
  double residue = -1;
};
