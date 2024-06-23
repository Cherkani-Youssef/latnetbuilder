#ifndef NETBUILDER__TASK__RANDOM_SEARCH_COLUMN_BY_COLUMN_H
#define NETBUILDER__TASK__RANDOM_SEARCH_COLUMN_BY_COLUMN_H

#include "netbuilder/Task/Search.h"
#include "latbuilder/LFSR258.h"
#include "netbuilder/GeneratingMatrix.h"
#include "netbuilder/Helpers/RankComputer.h"
#include <limits>

/**
 * This algorithm, described by Shin Harase in "Monte Carlo Methods and
 * Applications 22.4 (2016), pp. 349–357," efficiently searches for low WAFOM
 * point sets. It extends a point set \(P_d = \{\mathbf{x}_0, \dots,
 * \mathbf{x}_{2^d-1}\}\) for \(0 \leq d \leq k\), where each \(P_d\) is a
 * subset of \(P\) with \(|P|=2^k\). The algorithm evaluates numerous point sets
 * for each \(d\), retaining the best in terms of WAFOM.
 * 
 * This class constructs the columns of the generating matrix \(C\) column by
 * column randomly, ensuring the upper left (\(k \times k\)) submatrices are
 * invertible.
 */

namespace NetBuilder
{
    namespace Task
    {

        /** Class for Random Column By Column Search tasks.
         */
        template <NetConstruction NC, EmbeddingType ET, template <NetConstruction> class OBSERVER = MinimumObserver>
        class RandomSearchColumnByColumn : public Search<NC, ET, OBSERVER>
        {
            static_assert(NC == NetConstruction::EXPLICIT, "RandomSearchColumnByColumn only supports EXPLICIT method.");
            typedef NetConstructionTraits<NC> ConstructionMethod;
            using GenValue = typename ConstructionMethod::GenValue;
            using SizeParameter = typename ConstructionMethod::SizeParameter;

        public:
            /** Constructor.
             * @param dimension Dimension of the searched net.
             * @param sizeParameter Size parameter of the searched net.
             * @param nbTries Number of nets to evaluated.
             * @param figure Figure of merit used to compare nets.
             * @param verbose Verbosity level.
             * @param earlyAbortion Early-abortion switch. If true, the computations will be stopped if the net is worse than the best one so far.
             */
            RandomSearchColumnByColumn(Dimension dimension,
                                       SizeParameter sizeParameter,
                                       std::unique_ptr<FigureOfMerit::FigureOfMerit> figure,
                                       unsigned nbTries,
                                       int verbose = 0,
                                       bool earlyAbortion = false) : Search<NC, ET, OBSERVER>(dimension, sizeParameter, verbose, earlyAbortion),
                                                                     m_figure(std::move(figure)),
                                                                     m_nbTries(nbTries),
                                                                     m_sizeParameter(sizeParameter){};
            /** Constructor.
             * @param dimension Dimension of the searched net.
             * @param baseNet Net from which to start the search.
             * @param sizeParameter Size parameter of the net.
             * @param nbTries Number of nets to evaluated.
             * @param figure Figure of merit used to compare nets.
             * @param verbose Verbosity level.
             * @param earlyAbortion Early-abortion switch. If true, the computations will be stopped if the net is worse than the best one so far.
             */
            RandomSearchColumnByColumn(Dimension dimension,
                                       std::unique_ptr<DigitalNet<NC>> baseNet,
                                       SizeParameter sizeParameter,
                                       std::unique_ptr<FigureOfMerit::FigureOfMerit> figure,
                                       unsigned nbTries,
                                       int verbose = 0,
                                       bool earlyAbortion = false) : Search<NC, ET, OBSERVER>(dimension, std::move(baseNet), verbose, earlyAbortion),
                                                                     m_figure(std::move(figure)),
                                                                     m_nbTries(nbTries),
                                                                     m_sizeParameter(std::move(sizeParameter)){};

            /**
             * Default move constructor.
             * Deletes the implicit copy constructor.
             */
            RandomSearchColumnByColumn(RandomSearchColumnByColumn &&) = default;

            /**
             *  Returns information about the task
             */
            virtual std::string format() const override
            {
                std::string res;
                std::ostringstream stream;
                stream << Search<NC, ET, OBSERVER>::format();
                stream << "Exploration method: random Search Column By Column - " << m_nbTries << " samples" << std::endl;
                stream << "Figure of merit: " << m_figure->format() << std::endl;
                res += stream.str();
                stream.str(std::string());
                return res;
            }

            /**
             * Resets the search.
             */
            virtual void reset() override
            {
                Search<NC, ET, OBSERVER>::reset();
                this->m_figure->evaluator()->reset();
            }

            /**
             * Executes the search task.
             * The best net and merit value are set in the process.
             */
            virtual void execute() override
            {
                LatBuilder::LFSR258 rng;

                auto evaluator = this->m_figure->evaluator();
                // compute the merit of the base net is one was provided
                Real merit = 0;

                if (this->m_earlyAbortion)
                {
                    evaluator->onProgress().connect(boost::bind(&Search<NC, ET, OBSERVER>::Observer::onProgress, &this->observer(), boost::placeholders::_1));
                    evaluator->onAbort().connect(boost::bind(&Search<NC, ET, OBSERVER>::Observer::onAbort, &this->observer(), boost::placeholders::_1));
                }
                auto Oldnet = this->m_observer->bestNet();

                std::vector<GenValue> genVals, genVals_try;
                genVals.reserve(this->dimension());

                unsigned int startColIndex;
               
                // When we get a baseNet we need to store the generating matrices 
                if (this->nCols() != m_sizeParameter.second)
                {
                    for (Dimension dim = 0; dim < this->dimension(); ++dim)
                    {
                        GenValue genMat = Oldnet.generatingMatrix(dim);
                        genVals.push_back(std::move(genMat));
                    }
                    // genVals_try = genVals;
                    startColIndex = this->nCols();
                }
                else
                {
                    startColIndex = 0;
                }

                
                for (unsigned int colIndex = startColIndex + 1; colIndex <= m_sizeParameter.second; colIndex++)
                {

                    double lastmerit = std::numeric_limits<double>::infinity();

                    if (this->m_verbose > 0)
                    {
                        std::cout << "Column  " << colIndex << "/" << m_sizeParameter.second << std::endl;
                    }

                    for (unsigned int attempt = 1; attempt <= m_nbTries; ++attempt)
                    {
                        if (this->m_verbose > 0 && ((m_nbTries > 100 && attempt % 100 == 0) || (attempt % 10 == 0)))
                        {
                            std::cout << "Net " << attempt << "/" << m_nbTries << std::endl;
                        }

                        std::vector<GenValue> tmp_genVals = genVals;

                        RankComputer rankComputer(colIndex);
                        for (Dimension dim = 0; dim < this->dimension(); ++dim)
                        {

                            if (colIndex == 1)
                            {
                                GenValue col = GeneratingMatrix::createColumnMatrix(this->nRows(), colIndex, rng); // New column

                                GenValue subMat = col.upperLeftSubMatrix(colIndex, colIndex);

                                if (col(0, 0) == 0)
                                {
                                    subMat.flip(0, 0);
                                    col.flip(0, 0);
                                }

                                tmp_genVals.push_back(std::move(col));
                            }
                            else
                            {
                                GenValue col = GeneratingMatrix::createColumnMatrix(this->nRows(), colIndex, rng); // New column
                                tmp_genVals[dim].stackRight(std::move(col));

                                if (rankComputer.checkIfInvertible(tmp_genVals[dim].upperLeftSubMatrix(colIndex, colIndex)) == 0)
                                {
                                    tmp_genVals[dim].flip(colIndex - 1, colIndex - 1);
                                }
                            }
                        }

                        
                        SizeParameter sizee(this->nRows(), colIndex);
                        auto net = std::make_unique<DigitalNet<NC>>(this->m_dimension, sizee, tmp_genVals);

                        double merit = (*evaluator)(*net, this->m_verbose - 3);

                        if (merit < lastmerit)
                        {
                            lastmerit = merit;
                            genVals_try = std::move(tmp_genVals);
                            if (colIndex == m_sizeParameter.second)
                            {
                                this->m_observer->observe(std::move(net), merit);
                            }
                        }
                        if ((this->m_verbose > 0))
                        {

                            std::cout << "Net " << attempt << "/" << m_nbTries << " merit " << merit << std::endl;
                        }
                    }

                    // rankComputers = std::move(rankComputers_try);
                    // rng.jump();
                    genVals = std::move(genVals_try);
                }

                if (!this->m_observer->hasFoundNet())
                {
                    this->onFailedSearch()(*this);
                    return;
                }
                this->selectBestNet(this->m_observer->bestNet(), this->m_observer->bestMerit());

            }

            /**
             * {@inheritDoc}
             */
            virtual const FigureOfMerit::FigureOfMerit &figureOfMerit() const override
            {
                return *m_figure;
            }

        private:
            std::unique_ptr<FigureOfMerit::FigureOfMerit> m_figure;
            unsigned int m_nbTries;
            SizeParameter m_sizeParameter;
        };

    }
}

#endif
