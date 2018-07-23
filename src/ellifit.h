#pragma once

#include <vector>

#include <Eigen/Eigenvalues>

#include "ellipse.h"

inline Ellipse_t ellifit(const std::vector<Eigen::Vector2d>& points)
{
  /// http://www.sciencedirect.com/science/article/pii/S0031320312004694
  using namespace Eigen;
  Ellipse_t ellipse;
  const size_t N = points.size();
  if (N < 5)
  {
    return ellipse;
  }

  /// Step 1
  double mean_x = 0, mean_y = 0;
  double x_min = std::numeric_limits<int>::max(), y_min = std::numeric_limits<int>::max();
  double x_max = 0, y_max = 0;
  for (auto& point : points)
  {
    x_min = point.x() < x_min ? point.x() : x_min;
    x_max = point.x() > x_max ? point.x() : x_max;
    y_min = point.y() < y_min ? point.y() : y_min;
    y_max = point.y() > y_max ? point.y() : y_max;
  }

  mean_x = (x_min + x_max) / 2;
  mean_y = (y_min + y_max) / 2;

  /// Step 2
  MatrixXd X(N, 5);
  VectorXd Yb(N);

  // Construct matrices to solve for
  int i = 0;
  for (const auto& point : points)
  {
    const auto xb = point.x() - mean_x;
    const auto yb = point.y() - mean_y;
    X.row(i) << xb * xb, 2 * xb * yb, -2 * xb, -2 * yb, -1;
    Yb(i) = -yb * yb;
    i++;
  }

  /// Step 3

  /// Use QR-Decomposition instead of SigularValue-Decomposition, because it is faster
  /// and has better numerical stability
  auto qr = X.householderQr();
  const auto phi = qr.solve(Yb);

  /// Step 4
  const auto a = phi(0);
  const auto b = phi(1);
  const auto c = phi(2);
  const auto d = phi(3);
  const auto e = phi(4);

  const double xcyc_denominator = a - (b * b);
  const double x = (c - (d * b)) / xcyc_denominator;
  const double y = ((a * d) - (c * b)) / xcyc_denominator;

  const double ab_nominator = 2 * (e + (y * y) + (x * x * a) + (2 * b));
  const double first_denominator = (1 + a);
  const double second_denominator = std::sqrt((1 - a) * (1 - a) + (4 * b * b));

  ellipse.axes.x() = std::sqrt(ab_nominator / (first_denominator - second_denominator));
  ellipse.axes.y() = std::sqrt(ab_nominator / (first_denominator + second_denominator));

  ellipse.rotation = -0.5 * std::atan2(2 * b, 1 - a);

  ellipse.center.x() = x + mean_x;
  ellipse.center.y() = y + mean_y;

  ellipse.residue = (X * phi - Yb).norm() / Yb.norm();

  if (std::isnan(ellipse.axes.x())      //
      || std::isnan(ellipse.axes.y())   //
      || std::isnan(ellipse.center.x()) //
      || std::isnan(ellipse.center.y()) //
      || std::isnan(ellipse.rotation))
  {
    ellipse = Ellipse_t();
  }

  return ellipse;
}
