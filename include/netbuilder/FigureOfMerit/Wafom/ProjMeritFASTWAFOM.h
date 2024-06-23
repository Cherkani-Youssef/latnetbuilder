#ifndef NETBUILDER__FIGURE_OF_MERIT_BIT__PROJ_MERIT_FAST_WAFOM_H
#define NETBUILDER__FIGURE_OF_MERIT_BIT__PROJ_MERIT_FAST_WAFOM_H
#include <iostream>
#include <vector>
#include <numeric>
#include <chrono>
#include <future>
#include <thread>
#include <cstdint>
#include <bitset>

#include "netbuilder/FigureOfMerit/Wafom/GetColsReverseCache.h"
#include "netbuilder/FigureOfMerit/WeightedFigureOfMerit.h"
#include "netbuilder/Helpers/CBCCoordinateSet.h"


#include "latbuilder/Functor/LookUpTable.h"

/**
 * Template class representing a projection-dependent merit defined by the WAFOM of the projection.
 * 
 * WAFOM (Walsh Figure of Merit): WAFOM is a measure used in numerical integration to assess the 
 * quality of digital nets. It evaluates the uniformity of point distributions in a digital net by 
 * using Walsh functions. Lower WAFOM values indicate better uniformity and, consequently, better 
 * performance for numerical integration tasks.
 */

namespace NetBuilder
{
    namespace FigureOfMerit
    {

        class ProjMeritFASTWAFOM
        {
        public:
            typedef std::unique_ptr<LevelCombiner::LevelCombiner> pCombiner;

            /**
             * Constructs a projection-dependent merit based on the wafom of projections.
             * @param maxCardinal Maximum order of the subprojections.
             * @param combiner Not used. For sake of uniformity.
             */
            ProjMeritFASTWAFOM(uInteger w, uInteger q, const LookUpTable& table_c ,uInteger maxCardinal, pCombiner combiner = std::make_unique<LevelCombiner::LevelCombiner>()) : 
            m_w(w),
            m_q(q),
            m_maxCardinal(maxCardinal),
            m_table_c(table_c)
            {

                if (m_q <= 0) throw std::invalid_argument("q must be a positive integer that divides the number of rows and must not be zero");
                if (m_w % m_q != 0) throw std::invalid_argument("q must be a divisor of the number of rows");

            };

            /**
             * Returns the maximum order of the subprojections to take into account.
             */
            unsigned int maxCardinal() const { return m_maxCardinal; }

            /**
             * Output information about the figure of merit.
             */
            std::string format() const
            {
                std::string res;
                res += "Fast WAFOM based Figure of merit";
                res += "\nEmbedding type: Unilevel";
                return res;
            };
            Real operator()(const AbstractDigitalNet &net, const LatticeTester::Coordinates &projection) const
            {

                
                const uInteger k = net.numColumns();
                const uInteger w = net.numRows();
                const uInteger numPoints =  (1 << k);
                const uInteger l = w / m_q;
                double sum = 0.0;
       
                

                auto computeWAFOMRange = [&](int start, int end)
                {
                    GetColsReverseCache cache(net);
                    double localSum = 0.0;

                    uint64_t cachedCurPoint[projection.size()];
                    int bit;
                    double prod = 1.0;
                    int j = -1;

                    for (auto i = start; i < end; i++)
                    {
                       cache.getPointForProj( i, projection, cachedCurPoint);
                        prod = 1.0;
                        j = -1;
                        for (auto dim : projection)
                        {
                            j++;
                            uint64_t x_i = cachedCurPoint[j];

                            for (uInteger c = 1; c <= m_q + (w % m_q != 0 ? 1 : 0); ++c)
                            {
                                int d_c = cache.compute_d_c(x_i, c, w, m_q); // d_c is the integer representation of the c-th segment of the i-th row vector.
                                prod *= m_table_c.get( c , d_c ); 
                            }
                        }

                        localSum += prod  - 1.0;
                    }
                    return localSum;
                };

                 sum = computeWAFOMRange(0, numPoints);

                return (sum / numPoints);
            }

            
          

            CBCCoordinateSet projections(Dimension dimension) const
            {
                return CBCCoordinateSet(dimension, m_maxCardinal);
            }

        private:
            uInteger m_maxCardinal; // maximum order of subprojections to take into account
            uInteger m_w;
            uInteger m_q;
            const LookUpTable& m_table_c;
        };

    } // namespace FigureOfMerit

} // namespace NetBuilder

#endif
