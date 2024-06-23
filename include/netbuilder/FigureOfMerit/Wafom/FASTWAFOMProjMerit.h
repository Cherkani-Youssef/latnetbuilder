
#ifndef NETBUILDER__FIGURE_OF_MERIT_BIT__FAST_WAFOM_PROJ_MERIT_H
#define NETBUILDER__FIGURE_OF_MERIT_BIT__FAST_WAFOM_PROJ_MERIT_H
#include "netbuilder/FigureOfMerit/WeightedFigureOfMerit.h"
#include "netbuilder/FigureOfMerit/ProjectionDependentEvaluator.h"
#include "netbuilder/FigureOfMerit/LevelCombiner.h"
#include "netbuilder/FigureOfMerit/Wafom/GetColsReverseCache.h"
#include "latbuilder/Functor/LookUpTable.h"

#include <functional>
#include <stdexcept>
#include <vector>
#include <numeric>
#include <chrono>
#include <future>
#include <thread>
#include <cstdint>
#include <bitset>


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

        using LatticeTester::Coordinates;

        class FASTWAFOMProjMerit
        {

        public:
            typedef std::unique_ptr<LevelCombiner::LevelCombiner> pCombiner;

            /// Type of the merit value
            typedef double Merit;

            /// Type of the subprojection combination of merit
            typedef double SubProjCombination;
            /**
             * Constructs a projection-dependent merit based on the Wafom of projections.
             * @param maxCardinal Maximum order of the subprojections.
             * @param combiner Not used. For sake of uniformity.
             */
            FASTWAFOMProjMerit(uInteger w, uInteger q, const LookUpTable &table_c, uInteger maxCardinal, pCombiner combiner = std::make_unique<LevelCombiner::LevelCombiner>()) : m_w(w),
                                                                                                                                                                                  m_q(q),
                                                                                                                                                                                  m_maxCardinal(maxCardinal),
                                                                                                                                                                                  m_table_c(table_c)
            {

                if (m_q <= 0)
                    throw std::invalid_argument("q must be a positive integer that divides the number of rows and must not be zero");
                // if (m_w % m_q != 0) throw std::invalid_argument("q must be a divisor of the number of rows");
            };
            virtual ~FASTWAFOMProjMerit(){};

            unsigned int maxCardinal() const { return m_maxCardinal; }

            virtual std::string format() const
            {
                std::string res;
                res += "WAFOM based figure of merit";
                res += "\nEmbedding type: Unilevel";
                return res;
            };

            /**
             * Computes the projection-dependent merit of the net \c net for the given projection.
             * @param net Digital net to evaluate.
             * @param projection Projection to use.
             * @param maxMeritsSubProj Maximum of the wafom of the subprojections.
             */
            Real operator()(const AbstractDigitalNet &net, const LatticeTester::Coordinates &projection, SubProjCombination maxMeritsSubProj) const
            { // const uInteger k = net.numColumns();
                const uInteger k = net.generatingMatrix(0).nCols();
                const uInteger numIteration = (1 << k);
                const uInteger w = net.numRows();
                const uInteger numPoints = net.numPoints();
                const uInteger l = w / m_q;
                double sum = 0.0;

                // std::cout << numIteration <<" *********************** " << numPoints << std::endl;

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
                                prod *= m_table_c.get(c, d_c);               // lookup_tables[c-1][d_c];
                            }
                        }

                        localSum += prod  - 1.0 ;
                    }
                    return localSum;
                };

                sum = computeWAFOMRange(0, numIteration);

                return (sum / numPoints);
            }

            virtual Real combine(Merit merit, const AbstractDigitalNet &net, const LatticeTester::Coordinates &projection)
            {
                return (Real)merit;
            }

            /** Updates the combination of merit \c subProjCombination using \c merit.
             * @param merit Merit used to update.
             * @param subProjCombination  Combination of merit to update.
             */
            static void update(Merit merit, SubProjCombination &subProjCombination)
            {

                subProjCombination = merit; // std::max(merit, subProjCombination);///// ???????????
            }

            /** Initializes the combination of merit \c subProjCombination.
             * @param subProjCombination  Combination of merit to initialized.
             */
            static void setToZero(SubProjCombination &subProjCombination)
            {
                subProjCombination = 0;
            }

            /**
             * Determines the number of levels from a net.
             */
            static unsigned int numLevels(const AbstractDigitalNet &net)
            {
                return 1;
            }

            /**
             * Resize the combination of merit \c subProjCombination to match the number of levels \c numLevels.
             * @param subProjCombination  Combination of merit to resize.
             * @param numLevels Number of levels.
             */

            static void resize(SubProjCombination &subProjCombination, unsigned int numLevels = 1) {};

            /**
             * Returns the size of the combination of merit \c subProjCombination.
             */
            static unsigned int size(const SubProjCombination &subProjCombination)
            {
                return 1;
            }

        protected:
            unsigned int m_maxCardinal; // maximum order of subprojections to take into account
        private:
            uInteger m_w;
            uInteger m_q;
            const LookUpTable &m_table_c;
        };

        template <>
        class WeightedFigureOfMerit<FASTWAFOMProjMerit>::WeightedFigureOfMeritEvaluator : public ProjectionDependentEvaluator<FASTWAFOMProjMerit>
        {
        public:
            WeightedFigureOfMeritEvaluator(WeightedFigureOfMerit<FASTWAFOMProjMerit> *figure) : ProjectionDependentEvaluator(figure)
            {
            }
        };

    } //  namespace FigureOfMerit
} // namespace NetBuilder
#endif
