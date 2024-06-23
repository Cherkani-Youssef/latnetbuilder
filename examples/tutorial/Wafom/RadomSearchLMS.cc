
// #include "netbuilder/Task/WafomTasks/CoordinateColumnSearch.h"

#include "netbuilder/Task/WafomTasks/LMSFIX.h"

#include <iostream>
#include <memory>

#include "netbuilder/Types.h"
#include "netbuilder/DigitalNet.h"
#include "netbuilder/FigureOfMerit/WeightedFigureOfMerit.h"

#include "netbuilder/FigureOfMerit/Wafom/ProjMeritFASTWAFOM.h"
#include "netbuilder/FigureOfMerit/Wafom/FastWAFOMProjMerit.h"
#include "netbuilder/FigureOfMerit/Wafom/Wafom.h"
#include "netbuilder/FigureOfMerit/Wafom/FastWafom.h"

// #include "Path.h"

#include "netbuilder/Task/Eval.h"

#include "latbuilder/Util.h"

#include "latbuilder/Functor/LookUpTable.h"
// #include "matplotlibcpp.h"
#include "netbuilder/Types.h"
// #include "netbuilder/Helpers/JoeKuo.h"

#include "netbuilder/Task/RandomSearch.h"
#include "netbuilder/Task/FullCBCExplorer.h"
#include "netbuilder/Task/CBCSearch.h"
#include "netbuilder/Task/RandomCBCExplorer.h"

// #include "netbuilder/FigureOfMerit/Wafom/ProjMeritFASTWAFOM.h"
#include "netbuilder/FigureOfMerit/Wafom/FastWAFOMProjMerit.h"
#include "latticetester/ProductWeights.h"
#include "latticetester/OrderDependentWeights.h"

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
        w = std::atoi(argv[1]);
        k = std::atoi(argv[2]);
        q = std::atoi(argv[3]);
        dim = std::atoi(argv[4]);
        nbTries = std::atoi(argv[5]);
        h = std::atof(argv[6]);
        factor = std::atof(argv[7]); // 31 15 3 5 10 1.0 1.0
    }
    std::cout << "                                         w = " << w  << " q = " << q << " k = " << k << " dim " << dim << " nbTries " << nbTries << " h  " << h << " factor " << factor << std::endl;
    typedef typename NetConstructionTraits<NetConstruction::SOBOL>::SizeParameter SizeParameter;
    SizeParameter size(k);
    typedef typename NetConstructionTraits<NetConstruction::EXPLICIT>::SizeParameter SizeParameter2;
    SizeParameter2 size2(k, k);

    int l = w / q;
    LookUpTable table_c(w, q, l, h, factor);

    // auto figure = std::make_unique<FastWafom>(q, w, std::move(table_c));

    auto weights = std::make_unique<LatticeTester::ProductWeights>(.7);
    auto projDepMerit = std::make_unique<FASTWAFOMProjMerit>(w, q, std::move(table_c), 5);
    auto figure = std::make_unique<WeightedFigureOfMerit<FASTWAFOMProjMerit>>(std::numeric_limits<Real>::infinity(), std::move(weights), std::move(projDepMerit));

    //! [full-CBC_explorer]
    // auto task_0 = std::make_unique<RandomSearch<NetConstruction::SOBOL, EmbeddingType::UNILEVEL>>(dim, size, std::move(figure), nbTries);
    // auto explorer = std::make_unique<FullCBCExplorer<NetConstruction::SOBOL, EmbeddingType::UNILEVEL>>(dim, size);
    // auto task_0 = std::make_unique<CBCSearch<NetConstruction::SOBOL, EmbeddingType::UNILEVEL, FullCBCExplorer>>(dim, size, std::move(figure), std::move(explorer));

    //! [random-CBC_explorer]
    auto explorer = std::make_unique<RandomCBCExplorer<NetConstruction::SOBOL, EmbeddingType::UNILEVEL>>(dim, size, 70);
    auto task_0 = std::make_unique<CBCSearch<NetConstruction::SOBOL, EmbeddingType::UNILEVEL, RandomCBCExplorer>>(dim, size, std::move(figure), std::move(explorer));
    {
        auto start = std::chrono::high_resolution_clock::now();
        task_0->execute();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // std::cout << task->bestNet().format(OutputStyle::NET) << std::endl;
        std::cout << std::endl;
        std::cout << "----> task_0 Execution Time: " << duration.count() << " milliseconds" << " merit Value " << task_0->bestMeritValue() << std::endl;
    }
     std::cout  << " -------------------------------------------------------- " << std::endl;
    table_c = LookUpTable(w, q, l, h, factor);
    auto figure2 = std::make_unique<FastWafom>(q, w, std::move(table_c));
    auto baseNet = std::make_unique<DigitalNet<NetConstruction::SOBOL>>(task_0->bestNet());

    auto task = std::make_unique<RandomSearchLMS<NetConstruction::EXPLICIT, NetConstruction::SOBOL, EmbeddingType::UNILEVEL>>(dim, std::move(baseNet), size2, std::move(figure2), nbTries);
    
    //! [random-CBC_explorer]

    {
        auto start = std::chrono::high_resolution_clock::now();
        task->execute();
        //    std::cout  << " --------------- " << std::endl;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        // std::cout << task->bestNet().format(OutputStyle::NET) << std::endl;
        // std::cout << task->bestNet().format(OutputStyle::TERMINAL) << std::endl;
        std::cout << std::endl;
        std::cout << "Execution Time: " << duration.count() << " milliseconds" << " merit Value " << task->bestMeritValue() << std::endl;
    }
    // LookUpTable table_c2(w, q, l, h, factor);

    auto figure3 = std::make_unique<Wafom>(h, factor);
    auto evaluator = figure3->evaluator();
    auto start = std::chrono::high_resolution_clock::now();
    double value = (*evaluator)(task->bestNet());
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "  Wafom " << value << " log10 Wafom " << (std::log10(std::abs(value))) << " task duration " << duration.count() << std::endl;
    std::cout << " DIFF WAFOM " << (task->bestMeritValue() - value) << std::endl;
    std::cout << std::endl;
}
