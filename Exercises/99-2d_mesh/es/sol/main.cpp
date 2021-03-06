#include <iostream>
#include "src/mesh.hpp"

int main() {
  using namespace Geometry;

  std::string ans;

  Mesh m("mesh.msh");
  m.buildEdges();
  std::cout << "Print mesh" << std::endl << std::flush;
  std::cout << m << std::flush;

  std::cout << "\nWhich kind of BC are you interested in? " << std::flush;
  std::cin >> ans;
  m.bc( ans );
  
  return 0;
}
