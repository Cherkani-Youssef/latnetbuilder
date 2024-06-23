#include <iostream>
#include <memory>

#include "netbuilder/Types.h"
#include "netbuilder/DigitalNet.h"

#include "netbuilder/Task/Eval.h"

#include "latbuilder/Functor/LookUpTable.h"
// #include "matplotlibcpp.h"
#include "netbuilder/Types.h"

#include "netbuilder/Task/RandomSearch.h"

#include "netbuilder/FigureOfMerit/CoordUniformFigureOfMerit.h"
#include "latbuilder/Kernel/FastWAFOM.h"

#include "netbuilder/FigureOfMerit/Wafom/Wafom.h"
#include "netbuilder/FigureOfMerit/Wafom/FastWafom.h"

#include "latticetester/ProductWeights.h"
#include "latticetester/OrderDependentWeights.h"

// #include "matplotlibcpp.h"


using namespace NetBuilder;

using namespace NetBuilder::FigureOfMerit;
using namespace NetBuilder::Task;
// namespace plt = matplotlibcpp;

int main(int argc, char **argv)
{
    uInteger k, w, q, nbTries;
    Dimension dim;
    double h, factor;

    if (argc != 8)
    {
        // Parse command line arguments
        w = 15;
        k = w;
        q = 3;
        dim = 5;
        nbTries = 10;
        h = 1.0;
        factor = 2.0;
    }
    else
    {
        // Parse command line arguments
        k = std::atoi(argv[1]);
        w = k; // std::atoi(argv[2]);
        q = std::atoi(argv[3]);
        dim = std::atoi(argv[4]);
        nbTries = std::atoi(argv[5]);
        h = std::atof(argv[6]);
        factor = std::atof(argv[7]);

    }
    int l = w / q;

    typedef typename NetConstructionTraits<NetConstruction::SOBOL>::SizeParameter SizeParameter;
    SizeParameter size(w);
    LookUpTable table_c(w, q, l, h, factor);

   

    // auto weights = std::make_unique<LatticeTester::OrderDependentWeights>(.7);
    auto weights = std::make_unique<LatticeTester::ProductWeights>(.7);
    auto kernel = LatBuilder::Kernel::FastWAFOM(k, q, std::move(table_c));
    auto figure = std::make_unique<CoordUniformFigureOfMerit<LatBuilder::Kernel::FastWAFOM, EmbeddingType::UNILEVEL>>(std::move(weights), kernel);

    //! [RandomSearch]
    auto task = std::make_unique<RandomSearch<NetConstruction::SOBOL, EmbeddingType::UNILEVEL>>(dim, size, std::move(figure), nbTries);
    //! [RandomSearch]

    {
        auto start = std::chrono::high_resolution_clock::now();
        task->execute();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        std::cout << task->bestNet().format(OutputStyle::NET) << std::endl;
        // std::cout << task->bestNet().format(OutputStyle::TERMINAL) << std::endl;
        std::cout << std::endl;
        std::cout << "Execution Time: " << duration.count() << " milliseconds" << " merit Value " << task->bestMeritValue() << std::endl;
    }
    

    auto figure2 = std::make_unique<Wafom>(h, factor);
    auto evaluator = figure2->evaluator();
    auto start = std::chrono::high_resolution_clock::now();
    double value = (*evaluator)(task->bestNet());
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "  Wafom " << value << " log10 Wafom " << (std::log10(std::abs(value))) << " task duration " << duration.count() << std::endl;
    
}