
#ifndef LATBUILDER__KERNEL__WAFOM_H
#define LATBUILDER__KERNEL__WAFOM_H

#include "latbuilder/Kernel/FunctorAdaptor.h"
#include "latbuilder/Functor/WAFOM.h"


namespace LatBuilder { namespace Kernel {
class WAFOM : public FunctorAdaptor<Functor::WAFOM> {
public:
   WAFOM(unsigned int outdigits,  double h, double factor):
          FunctorAdaptor<Functor>(Functor(outdigits, h, factor))
   {}
   static constexpr Real CUPower = 2;

};
}
} 

#endif
