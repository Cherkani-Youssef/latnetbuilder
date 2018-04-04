// This file is part of Lattice Builder.
//
// Copyright (C) 2012-2016  Pierre L'Ecuyer and Universite de Montreal
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "latbuilder/Storage.h"
#include "latbuilder/SizeParam.h"

#include <iostream>

using namespace LatBuilder;

//! [all]
template <LatticeType LA, LatEmbed LAT, Compress COMP>
void test(typename LatticeTraits<LA>::Modulus modulus)
{
   SizeParam<LA, LAT> size(modulus);
   Storage<LA, LAT, COMP> storage(size);
   std::cout << "storage name: " << storage.name() << std::endl;
   std::cout << "  size parameter: " << storage.sizeParam() << std::endl;
   std::cout << "  virtual size:   " << storage.virtualSize() << std::endl;
   std::cout << "  actual size:    " << storage.size() << std::endl;
}

int main()
{

   uInteger n = 16;
   Polynomial P = PolynomialFromInt(7);
   test<LatticeType::ORDINARY, LatEmbed::SIMPLE, Compress::NONE>(n);
   test<LatticeType::ORDINARY, LatEmbed::EMBEDDED, Compress::NONE>(n);
   test<LatticeType::ORDINARY, LatEmbed::SIMPLE, Compress::SYMMETRIC>(n);
   test<LatticeType::ORDINARY, LatEmbed::EMBEDDED, Compress::SYMMETRIC>(n);

   test<LatticeType::POLYNOMIAL, LatEmbed::SIMPLE, Compress::NONE>(P);
   test<LatticeType::POLYNOMIAL, LatEmbed::EMBEDDED, Compress::NONE>(P);

   return 0;
}
//! [all]
