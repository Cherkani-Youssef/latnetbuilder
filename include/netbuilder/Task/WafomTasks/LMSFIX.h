#ifndef NETBUILDER__TASK__RANDOM_SEARCH_LMS_H
#define NETBUILDER__TASK__RANDOM_SEARCH_LMS_H

#include "netbuilder/Task/Search.h"
#include "latbuilder/LFSR258.h"

namespace NetBuilder
{
    namespace Task
    {

        template <NetConstruction NC, NetConstruction NCBaseNet, EmbeddingType ET, template <NetConstruction> class OBSERVER = MinimumObserver>
        class RandomSearchLMS : public Search<NC, ET, OBSERVER>
        {
            LatBuilder::LFSR258 rng;
            typedef NetConstructionTraits<NC> ConstructionMethod;
            // typedef NetConstructionTraits<NC> ConstructionMethod;
            using GenValue = typename ConstructionMethod::GenValue;
            using SizeParameter = typename ConstructionMethod::SizeParameter;

        public:
            /** Constructor.
             * @param dimension Dimension of the searched net.
             * @param baseNet Net from which to start the search.
             * @param sizeParameter Size parameter of the net.
             * @param nbTries Number of nets to evaluated.
             * @param figure Figure of merit used to compare nets.
             * @param verbose Verbosity level.
             * @param earlyAbortion Early-abortion switch. If true, the computations will be stopped if the net is worse than the best one so far.
             */

            RandomSearchLMS(Dimension dimension,
                            std::unique_ptr<DigitalNet<NCBaseNet>> baseNet,
                            typename NetConstructionTraits<NC>::SizeParameter sizeParameter,
                            std::unique_ptr<FigureOfMerit::FigureOfMerit> figure,
                            unsigned nbTries,
                            int verbose = 0,
                            bool earlyAbortion = false) : Search<NC, ET, OBSERVER>(dimension, sizeParameter, verbose, earlyAbortion),
                                                          m_figure(std::move(figure)),
                                                          m_nbTries(nbTries),
                                                          m_sizeParameter(std::move(sizeParameter)),
                                                          m_randomGenValueGenerator(this->m_sizeParameter),
                                                          m_baseNet(std::move(baseNet)){};

            /**
             * Default move constructor.
             * Deletes the implicit copy constructor.
             */
            RandomSearchLMS(RandomSearchLMS &&) = default;

            /**
             *  Returns information about the task
             */
            virtual std::string format() const override
            {
                std::string res;
                std::ostringstream stream;
                stream << Search<NC, ET, OBSERVER>::format();
                stream << "Exploration method: random - " << m_nbTries << " samples" << std::endl;
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

                unsigned int w = this->nRows();
                unsigned int k = this->nCols();

                auto evaluator = this->m_figure->evaluator();
                // auto baseNet = this->m_observer->bestNet();

                if (this->m_earlyAbortion)
                {
                    evaluator->onProgress().connect(boost::bind(&Search<NC, ET, OBSERVER>::Observer::onProgress, &this->observer(), boost::placeholders::_1));
                    evaluator->onAbort().connect(boost::bind(&Search<NC, ET, OBSERVER>::Observer::onAbort, &this->observer(), boost::placeholders::_1));
                }

                for (unsigned int attempt = 1; attempt <= m_nbTries; ++attempt)
                {

                    if (this->m_verbose > 0 && ((m_nbTries > 100 && attempt % 100 == 0) || (attempt % 10 == 0)))
                    {
                        std::cout << "Net " << attempt << "/" << m_nbTries << std::endl;
                    }

                    
                     std::vector<typename ConstructionMethod::GenValue> genVals;
                     genVals.reserve(this->dimension());
                    for (Dimension dim = 0; dim < m_baseNet->dimension(); ++dim)
                    {
                        auto C = m_baseNet->generatingMatrix(dim);
                        auto L = (GeneratingMatrix::createRandomLowerTriangularMatrix(w, w, rng));
                        auto C_tilde = (L * C); 
                        genVals.push_back(std::move(C_tilde));
                    }
                   
                    auto net = std::make_unique<DigitalNet<NC>>(this->m_dimension, this->m_sizeParameter, std::move(genVals));

                    double merit = (*evaluator)(*net, this->m_verbose - 3);

                    this->m_observer->observe(std::move(net), merit);
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
            typename ConstructionMethod::template RandomGenValueGenerator<ET> m_randomGenValueGenerator;

            std::unique_ptr<DigitalNet<NCBaseNet>> m_baseNet;
        };

    }
}

#endif
