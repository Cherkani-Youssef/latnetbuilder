

#include <iostream>
#include <memory>

#include "netbuilder/Types.h"
#include "netbuilder/DigitalNet.h"
#include "netbuilder/FigureOfMerit/WeightedFigureOfMerit.h"
#include "netbuilder/FigureOfMerit/Wafom/Wafom.h"
#include "netbuilder/FigureOfMerit/TValue.h"
#include "netbuilder/FigureOfMerit/Wafom/FastWafom.h"

#include "netbuilder/Task/WafomTasks/RandomSearchColumnByColumn.h"

#include "latticetester/ProductWeights.h"

// #include "Path.h"
// #include "WafomHelpers.h"
#include "netbuilder/Task/Eval.h"

#include "latbuilder/Util.h"

#include "latbuilder/Functor/LookUpTable.h"

using namespace NetBuilder;

using namespace NetBuilder::FigureOfMerit;
using namespace NetBuilder::Task;
int main(int argc, char **argv)
{

    uInteger k, w, q, nbTries;
    Dimension dim;
    double h, factor;

    if (argc != 8)
    {
        // Parse command line arguments
        w = 31;
        k = 15;
        q = 3;
        dim = 5;
        nbTries = 10;

        h = 1.0;
        factor = 2.0;
    }
    else
    {
        // Parse command line arguments
        w = std::atoi(argv[1]);
        k = std::atoi(argv[2]);
        q = std::atoi(argv[3]);
        dim = std::atoi(argv[4]);
        nbTries = std::atoi(argv[5]);
        h = std::atof(argv[6]);
        factor = std::atof(argv[7]);
    }

    // typedef typename NetConstructionTraits<NetConstruction::SOBOL>::SizeParameter SizeParameter;
    // SizeParameter size(k);

    typedef typename NetConstructionTraits<NetConstruction::EXPLICIT>::SizeParameter SizeParameter;

    SizeParameter size(w, k);

    int l = w / q;
    LookUpTable table_c(w, q, l, h, factor);

    auto figure = std::make_unique<FastWafom>(q, w, std::move(table_c));

    //! [random_search_task]
    auto task = std::make_unique<RandomSearchColumnByColumn<NetConstruction::EXPLICIT, EmbeddingType::UNILEVEL>>(dim, size, std::move(figure), nbTries);
    // auto start = std::chrono::high_resolution_clock::now();
    // task->execute();
    // auto end = std::chrono::high_resolution_clock::now();

    // std::cout << task->bestNet().format(OutputStyle::NET) << std::endl;

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
    // LookUpTable table_c2(w, q, l, h, factor);

    auto figure2 = std::make_unique<Wafom>(h, factor);
    auto evaluator = figure2->evaluator();
    auto start = std::chrono::high_resolution_clock::now();
    double value = (*evaluator)(task->bestNet());
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "  Wafom " << value << " log10 Wafom " << (std::log10(std::abs(value))) << " task duration " << duration.count() << std::endl;
    std::cout << " DIFF WAFOM" << (task->bestMeritValue() - value) << std::endl;
    std::cout << std::endl;
}