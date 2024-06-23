
#ifndef LATBUILDER__KERNEL__FastWAFOM_H
#define LATBUILDER__KERNEL__FastWAFOM_H

#include "latbuilder/Kernel/FunctorAdaptor.h"
#include "latbuilder/Functor/FastWAFOM.h"

namespace LatBuilder
{
   namespace Kernel
   {
      class FastWAFOM : public FunctorAdaptor<Functor::FastWAFOM>
      {
      public:
         FastWAFOM(uInteger outdigits, uInteger q, const LookUpTable &table_c) : FunctorAdaptor<Functor>(Functor(outdigits, q, table_c))
         {
         }

         static constexpr Real CUPower = 1;
      };
   }
}

#endif
