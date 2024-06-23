#ifndef NET_BUILDER__FIGURE_OF_MERIT__Wafom_H
#define NET_BUILDER__FIGURE_OF_MERIT__Wafom_H

#include <iostream>
#include <vector>
#include <numeric>
#include <chrono>
#include <future>
#include <thread>
#include <cstdint>
#include <bitset>

#include <mutex>
#include <cmath>

#include "netbuilder/FigureOfMerit/FigureOfMerit.h"
#include "netbuilder/FigureOfMerit/Wafom/GetColsReverseCache.h"

/**
 * This class calculates the Walsh Figure of Merit (WAFOM) for a Digital Net Base 2.
 * WAFOM is a discrepancy measure used for evaluating point sets. It is defined as:
 *
 * WAFOM(P) = -1 + (1 / |P|) * Σ_{B ∈ P} {Π_{i=1}^{s} Π_{j=1}^{w} [(1 + (-1)^{b_{i,j}} * 2^{-j}) ]}
 *
 * Where |P| is the number of points in P, s is the dimension of the points, w is the precision,
 * and b_{i,j} is the  l-th bit of the j-th coordinate of a point x \in P.
 *
 * This discrepancy measure, specific to digital nets, is designed to provide excellent results for alpha-smooth functions
 * (refer to J. Dick. “On Quasi-Monte Carlo Rules Achieving Higher Order Convergence”. In: Monte Carlo and Quasi-Monte Carlo Methods 2008).
 */

namespace NetBuilder
{
    namespace FigureOfMerit
    {

        class Wafom : public FigureOfMerit
        {

        public:
            /**
             *
             *
             * @param digitalNet The Digital Net Base 2 instance.
             * @param outDigits  The precision of the points.
             * @param h          The parameter defining Matsumoto's or Yoshiki's definition.
             * @param factor     set to to 1 for wafom set to 2 to get Wafom for RMSE
             */
            Wafom(double h, double factor) : h(h), factor(factor)
            {
                if (h != 0 && h != 1)
                {
                    throw std::invalid_argument("h must be either 0 or 1: for Matsumoto's definition set to 0, and for Yoshiki's set to 1");
                }
                if (factor != 1 && factor != 2)
                {
                    throw std::invalid_argument("factor must be set to 1 for wafom, and set to 2 for Wafom for RMSE");
                }
            }
            std::unique_ptr<FigureOfMeritEvaluator> createEvaluator() override
            {
                return std::make_unique<WafomEvaluator>(this, h, factor);
            }

            virtual std::string format() const override
            {
                std::string res;
                res += "Wafom based figure of merit";
                res += "\nEmbedding type: Unilevel";
                return res;
            }

        private:
            double h;
            double factor;

            /**
             * Evaluator for the Walsh figure of merit.
             */
            class WafomEvaluator : public FigureOfMeritEvaluator
            {
            public:
                /**
                 * Constructor.
                 * @param figure Pointer to the figure of merit.
                 */
                WafomEvaluator(Wafom *figure, double h, double factor) : m_dimension(0),
                                                                         m_figure(figure), h(h), factor(factor){

                                                                                                 };

                /**
                 * Computes the figure of merit for the given \c net for all the dimensions (full computation).
                 * @param net Net to evaluate.
                 * @param verbose Verbosity level.
                 */
                virtual MeritValue operator()(const AbstractDigitalNet &net, int verbose = 0) override
                {

                    double sum = 0;
                    const int w = net.numRows();
                    const int numPoints = net.numPoints();
                    const int k = net.numColumns();
                    const int s = net.dimension();
                    GetColsReverseCache cache(net);

                    // Define a lambda function to compute WAFOM for a range of points
                    auto computeWAFOMRange = [&](int start, int end)
                    {
                        double localSum = 0.0;
                        double prod = 1.0;
                        size_t dim = net.dimension();
                        uint64_t cachedCurPoint[dim];
                        int bit;

                        for (int i = start; i < end; ++i)
                        {
                            prod = 1.0;

                            cache.getPoint(i, s, cachedCurPoint);
                            for (int j = 0; j < s; ++j)
                            {

                                for (int l = 1; l <= w; ++l)
                                {
                                    bit = ((cachedCurPoint[j] >> (w - l)) & 1);

                                    double two_exponent = std::pow(2.0, -(factor * (l + h)));
                                    prod *= (1 + (1 - 2 * bit) * two_exponent);
                                }
                            }
                            localSum += prod - 1.0;
                        }
                        return localSum;
                    };

                    sum = computeWAFOMRange(0, numPoints);

                    return (sum / numPoints);
                }

                /**
                 * Resets the evaluator and prepare it to evaluate a new net.
                 */
                virtual void reset() override
                {
                    return;
                }

            private:
                Dimension m_dimension;
                Wafom *m_figure;
                double h;
                double factor;
            };
        };

    }
}

#endif
