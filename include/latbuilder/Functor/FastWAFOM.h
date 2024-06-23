#ifndef LATBUILDER__FUNCTOR__FastWAFOM_H
#define LATBUILDER__FUNCTOR__FastWAFOM_H

#include <iostream>
#include <sstream>
#include <vector>
#include <cstdint>
#include "latbuilder/Types.h"
#include "latbuilder/Functor/LookUpTable.h"
// #include "netbuilder/FigureOfMerit/Wafom/GetColsReverseCache.h"

/**
 * One-dimensional merit function for the $Wafom$ discrepancy.
 *
 * This merit function is defined as
 * \f[
 *    \omega(x) =
 *    $\prod_{1\leq c \leq q} table)c[d(c_i)]$
 * \f]
 *  A coordinate of a point is expressed as \(X^i = d_1^i, \ldots, d_q^i \)
 *  where each segment is of length less or equal to l bits
 */

namespace LatBuilder
{
    namespace Functor
    {

        class FastWAFOM
        {
        public:
            // Assuming Real and uInteger are defined in "latbuilder/Types.h"
            typedef Real value_type;
            typedef Real result_type;

            // Pass LookupTable by reference to avoid unnecessary copies
            FastWAFOM(uInteger outdigits, uInteger q, const LookUpTable &table_c) : m_outdigits(outdigits),
                                                                                    m_q(q),
                                                                                    m_table_c(std::move(table_c))
            {

                if (m_q <= 0)
                    throw std::invalid_argument("q must be a positive integer that divides the number of rows and must not be zero");
                // if (m_outdigits % m_q != 0) throw std::invalid_argument("q must be a divisor of the number of rows");
            }

            bool symmetric() const { return false; }

            static constexpr Compress suggestedCompression() { return Compress::NONE; }

            // Correct operator() with const qualifier and static_cast
            result_type operator()(const value_type &x, uInteger n = 0) const
            {
                result_type product = 1.0;

                double exponent = static_cast<double>(1 << m_outdigits);
                uint64_t intRepresentation = static_cast<uint64_t>((x * exponent));
                for (uInteger c = 1; c <= m_q + (m_outdigits % m_q != 0 ? 1 : 0); ++c)
                {
                    uInteger d_c = compute_d_c(intRepresentation, c, m_outdigits, m_q);
                    product *= m_table_c.get(c, d_c);
                }
                return product;
            }
            /**
             * @param uint64_t x_i integer representation of the point,
             * @param c is the  index of the segment
             * @param w is the total number of bits
             * @param q number of segments
             */
            uInteger compute_d_c(uint64_t x_i, uInteger c, uInteger w, uInteger q) const
            {
                uInteger l = w / q;                              // Length of each full segment
                uInteger remainder = w % q;                      // Length of the remainder segment, if any
                uInteger segments = q + (remainder > 0 ? 1 : 0); // Total number of segments including remainder if exists

                // Boundary check for c
                if (c < 1 || c > segments)
                {
                    throw std::invalid_argument("Segment index c is out of bounds.");
                }

                if (c <= q)
                {
                    uInteger startBitFromMsb = (c * l) - 1;
                    uInteger shiftAmount = w - startBitFromMsb - 1;
                    return (x_i >> shiftAmount) & ((1 << l) - 1);
                }
                else
                {
                    // Handle the remainder segment correctly
                    return x_i & ((1 << remainder) - 1); // Extract the remainder segment value
                }
            }

            std::string name() const { return "FastWAFOM"; }

        private:
            uInteger m_outdigits;
            uInteger m_q;
            const LookUpTable &m_table_c;
        };

        inline std::ostream &operator<<(std::ostream &os, const FastWAFOM &functor)
        {
            return os << functor.name();
        }

    }
}

#endif // LATBUILDER__FUNCTOR__FastWAFOM_H
