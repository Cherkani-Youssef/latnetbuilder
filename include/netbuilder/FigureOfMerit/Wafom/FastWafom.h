
#ifndef NET_BUILDER__FIGURE_OF_MERIT__FASTWAFOM_H
#define NET_BUILDER__FIGURE_OF_MERIT__FASTWAFOM_H
#include <iostream>
#include <vector>
#include <array>
#include <numeric>
#include <chrono>
#include <future>
#include <thread>
#include <cstdint>
#include <bitset>

#include "netbuilder/FigureOfMerit/FigureOfMerit.h"

#include "netbuilder/FigureOfMerit/Wafom/GetColsReverseCache.h"
#include "latbuilder/Storage.h"
#include "latbuilder/Functor/LookUpTable.h"

/**
 * This class computes the Walsh Figure of Merit (WAFOM) for a Digital Net in
 * Base 2 using a lookup table to accelerate the computation.
 *
 * The WAFOM is a measure used in numerical integration to assess the 
 * quality of digital nets. It evaluates the uniformity of point distributions in a digital net by 
 * using Walsh functions. Lower WAFOM values indicate better uniformity and, consequently, better 
 * performance for numerical integration tasks.
 *
 * For Matsumoto's definition
 *
 * WAFOM(P) = -1 + \frac{1}{|P|} \left( \sum_{\mathbf{x}_i \in P}
 * \prod_{j=1}^{s} \prod_{l=1}^{w} \left( 1 + \eta(x_{i,j,l}) 2^{-l} \right)
 * \right)
 *
 * The same formula applies to Yoshiki's definition; simply replace \( l \) with
 * \( l + 1 \).
 *
 * WAFOM(P) = -1 + \frac{1}{|P|} \left( \sum_{\mathbf{x}_i \in P}
 * \prod_{j=1}^{s} \prod_{l=1}^{w} \left( 1 + \eta(x_{i,j,l}) 2^{-(l+1)} \right)
 * \right) .
 *
 * For Goda's definition:
 *
 * \[ \mathcal{W}(P, \mu) = \sqrt{-1 + \frac{1}{|P|} \left( \sum_{\mathbf{x}_i
 * \in P} \prod_{j=1}^{s} \prod_{l=1}^{w} \left( 1 + \eta(x_{i,j,l}) 2^{-2l}
 * \right) \right)} \] The same formula applies to Yoshiki's definition; simply
 * replace \( l \) with \( l + 1 \).
 *
 * The same formula applies to Yoshiki's definition; simply replace \( l \) with
 * \( l + 1 \).
 *
 * \[ \mathcal{W}(P, \mu + h) = \sqrt{-1 + \frac{1}{|P|} \left(
 * \sum_{\mathbf{x}_i \in P} \prod_{j=1}^{s} \prod_{l=1}^{w} \left( 1 +
 * \eta(x_{i,j,l}) 2^{-2(l+1)} \right) \right)} \]
 *
 * The method, proposed by Shin Harase in "A search for extensible low-WAFOM
 * point sets," Monte Carlo Methods and Applications 22.4 (2016), pp. 349–357,
 * involves dividing a coordinate of a point with 'w' precision into 'q'
 * subparts, each containing 'l' elements if q is a divisor of the number of
 * 'w'. When 'q' is not a divisor of the number of output digits then will
 * another be a segment of lehgth t where \(t = w \mod q\) .
 *
 * A coordinate of a point is expressed as \(X^i = d_1^i, \ldots, d_q^i,
 * d_{q+1}^i\), where for \(1 \leq c \leq q\), the length of each segment is \(l
 * = w \div q\), and the remainder \(d_{q+1}\) is of length \(t = w \mod q\).
 * Each segment \(d_c^i\) is represented as \(x_{i, (c-1)*l+1}, \ldots, x_{i,
 * c*l}\).
 *
 * Instead of computing \(\prod_{j=1}^{w} [1 + (-1)^{(-1)^{x_{i,j}}} 2^{-j}] -
 * 1\), we use \(\prod_{1 \leq c \leq q} \text{table}_c[d_c^i]\), where
 * \(\text{table}_c[d_c^i]\) are precomputed values: \(\text{table}_c[d_c^i] =
 * \prod_{1 \leq j \leq l} (1 + (-1)^{x_{i, (c-1)*l + j}} \cdot 2^{-((c-1)*l + j
 * + 1)})\). This significantly reduces computation time.
 *
 * We use Kahan summation algorithm improves numerical accuracy by compensating
 * for floating-point errors, reducing round-off error accumulation.
 **/

namespace NetBuilder
{
    namespace FigureOfMerit
    {

        class FastWafom : public FigureOfMerit
        {
        public:
            FastWafom(uInteger q, uInteger w, const LookUpTable &table_c) : m_q(std::move(q)), m_numRows(std::move(w)), m_table_c(table_c)
            {
                if (m_q <= 0)
                    throw std::invalid_argument("q must be a positive integer that divides the number of rows and must not be zero ");
                // if (m_numRows % m_q != 0)
                // {
                //     throw std::invalid_argument("q must be a divisor of the number of rows 0000");
                // }
            }

            std::unique_ptr<FigureOfMeritEvaluator> createEvaluator() override
            {
                return std::make_unique<FastWafomEvaluator>(this, m_q, m_numRows, m_table_c);
            }

            virtual std::string format() const override
            {

                std::string res;
                res += "Fast Wafom based figure of merit";
                res += "\nEmbedding type: Unilevel";
                return res;
            }

        private:
            class FastWafomEvaluator : public FigureOfMeritEvaluator
            {
            public:
                FastWafomEvaluator(FastWafom *figure, uInteger q, uInteger w, const LookUpTable &table_c)
                    : m_dimension(0), m_figure(figure), m_q(q), m_numRows(w), m_table_c(table_c) {}

                virtual MeritValue operator()(const AbstractDigitalNet &net, int verbose = 0) override
                {
                    const uInteger numPoints = net.numPoints();
                    const uInteger w = net.numRows();
                    const uInteger dim = net.dimension();
                    const uInteger l = w / m_q;
                    double sum = 0.0;

                    // Initialize cache with net
                    GetColsReverseCache cache(net);

                    uint64_t point[dim];

                    auto computeWAFOMRange = [&](uInteger start, uInteger end)
                    {
                        double localSum = 0.0;
                        double product = 1.0;

                        // Iterate over points to compute the WAFOM value
                        for (uInteger i = start; i < end; ++i)
                        {
                            product = 1.0;
                            cache.getPoint(i, dim, point);

                            for (uInteger j = 0; j < dim; ++j)
                            {
                                uint64_t u_i_j = point[j];
                                for (uInteger c = 1; c <= m_q + (w % m_q != 0 ? 1 : 0); ++c)
                                {
                                    uInteger d_c = cache.compute_d_c(u_i_j, c, w, m_q); 
                                    product *= m_table_c.get(c, d_c);
                                }
                            }
                            localSum += product - 1.0;
                        }
                        return localSum;
                    };

                    sum = computeWAFOMRange(0, numPoints);

                    return sum / numPoints;
                }

                virtual void reset() override {}

            private:
                Dimension m_dimension;
                FastWafom *m_figure;
                uInteger m_q;
                uInteger m_numRows;
                const LookUpTable &m_table_c;
            };

            uInteger m_q;
            uInteger m_numRows;
            const LookUpTable &m_table_c;
        };

    } // namespace FigureOfMerit
} // namespace NetBuilder

#endif
