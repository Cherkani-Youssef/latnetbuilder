#ifndef LATBUILDER__FUNCTOR__WAFOM_H
#define LATBUILDER__FUNCTOR__WAFOM_H

#include <iostream>

#include "latbuilder/Types.h"
#include <sstream>
#include <vector>

/**
 * One-dimensional merit function for the $Wafom$ discrepancy.
 *
 * This merit function is defined as
 * \f[
 *    \omega(x) = \prod_{l=1}^{w} \left( 1 + \eta(x_{i,j,l}) 2^{-(l+1)}
 * \f]
 *   where $\eta =  (−1)^{x_i}$
 */

namespace LatBuilder
{
   namespace Functor
   {
      class WAFOM
      {
      public:
         typedef Real value_type;
         typedef Real result_type;

         /**
          * Constructor.
          */
         WAFOM(unsigned int outdigits, double h, double factor) : m_outdigits(outdigits)
         {
         }

         bool symmetric() const
         {
            return false;
         }

         static constexpr Compress suggestedCompression()
         {
            return Compress::NONE;
         }

         result_type operator()(const value_type &x, unsigned long outdigits) const
         {
            double exponent = (1 << m_outdigits);
            int intRepresenatation = static_cast<int>(std::round(x * exponent));
            double prod = 1.0;
            int bit;
            for (int l = 1; l <= m_outdigits; ++l)
            { // Iterating from MSB to LSB
               bit = ((intRepresenatation >> (m_outdigits - l)) & 1);
               double two_exponent = std::pow(2.0, -(factor * (l + h)));
               prod *= (1 + (1 - 2 * bit) * two_exponent);
            }
            return prod;
         }

         std::string name() const
         {
            return " WAFOM ";
         }

      private:
         unsigned int m_outdigits;
         double h;
         double factor;
      };

      inline std::ostream &operator<<(std::ostream &os, const WAFOM &functor)
      {
         return os << functor.name();
      }

   }
   static constexpr Real CUPower = 2;
}

#endif
