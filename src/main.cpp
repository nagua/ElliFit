#include <iostream>
#include <random>
#include <string>

#include "ellifit.h"

int main(int argc, char* argv[])
{
  std::vector<Eigen::Vector2d> points;
  const int numberOfArguments = argc - 1;

  if (numberOfArguments >= 10 && numberOfArguments % 2 == 0)
  {
    std::cout << "Using user supplied values for ellipse fitting" << std::endl;
    points.resize(numberOfArguments / 2);
    for (int i = 0; i < numberOfArguments; i += 2)
    {
      auto& point = points[i / 2];
      point.x() = std::stod(argv[i + 1]);
      point.y() = std::stod(argv[i + 2]);
    }
  }
  else
  {
    std::cout << "Generating 20 random points to do ellipse fitting" << std::endl;
    points.resize(20);
    std::random_device rd;
    std::mt19937 gen(rd());
    const std::uniform_real_distribution<double> distribution(-20, 20);

    for (auto& point : points)
    {
      point.x() = distribution(gen);
      point.y() = distribution(gen);
    }
  }

  const auto ellipse = ellifit(points);


  if (ellipse.residue == -1)
  {
    std::cout << "Could not fit a ellipse" << std::endl;
  }
  else
  {
    std::cout << "Center:\t\t " << ellipse.center.transpose() << std::endl;
    std::cout << "Axes:\t\t " << ellipse.axes.transpose() << std::endl;
    std::cout << "Rotation:\t " << ellipse.rotation << std::endl;
    std::cout << "Residue:\t " << ellipse.residue << std::endl;
  }

  return 0;
}
