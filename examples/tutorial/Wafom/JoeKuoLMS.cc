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

#include "netbuilder/Task/Eval.h"

#include "latbuilder/Util.h"

#include "latbuilder/Functor/LookUpTable.h"

#include "netbuilder/Types.h"

#include "netbuilder/Helpers/JoeKuo.h"

#include "netbuilder/Task/RandomSearch.h"
#include "netbuilder/Task/FullCBCExplorer.h"
#include "netbuilder/Task/CBCSearch.h"
#include "netbuilder/Task/RandomCBCExplorer.h"


#include "netbuilder/FigureOfMerit/Wafom/FastWAFOMProjMerit.h"
#include "latticetester/ProductWeights.h"
#include "latticetester/OrderDependentWeights.h"

using namespace NetBuilder;

using namespace NetBuilder::FigureOfMerit;
using namespace NetBuilder::Task;
// namespace plt = matplotlibcpp;
#include "Path.h"
using JoeKuo::createJoeKuoSobolNet;

int main(int argc, char **argv)
{
    SET_PATH_TO_LATNETBUILDER_FOR_EXAMPLES();

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
        factor = std::atof(argv[7]); // 15 15 3 5 10 1.0 1.0
    }
    int l = w / q;
    LookUpTable table_c(w, q, l, h, factor);

    typedef typename NetConstructionTraits<NetConstruction::EXPLICIT>::SizeParameter SizeParameter;
    SizeParameter size(w,k);
std::cout << "                                         w = " << w  << " q = " << q << " k = " << k << " dim " << dim << " nbTries " << nbTries << " h  " << h << " factor " << factor << std::endl;
// auto baseNet = std::make_unique<DigitalNet<NetConstruction::SOBOL>>(task_0->bestNet());

    auto baseNet = std::make_unique<DigitalNet<NetConstruction::SOBOL>>(createJoeKuoSobolNet(dim, k));
    // std::cout << baseNet->format(OutputStyle::NET) << std::endl;
    auto figure2 = std::make_unique<FastWafom>(q, w, std::move(table_c));

   { auto fig = std::make_unique<FastWafom>(q, w, std::move(table_c));

      auto evaluator = fig->evaluator();
    double value = (*evaluator)(*baseNet);
    std::cout << "  Wafom " << value << " log10 Wafom " << (std::log10(std::abs(value)))<< std::endl;
   }
    auto task = std::make_unique<RandomSearchLMS<NetConstruction::EXPLICIT, NetConstruction::SOBOL, EmbeddingType::UNILEVEL>>(dim, std::move(baseNet), size, std::move(figure2), nbTries);
    //! [random-CBC_explorer]

    {
        auto start = std::chrono::high_resolution_clock::now();
        task->execute();
        //    std::cout  << " --------------- " << std::endl;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << " ------------------------------------------------------------------------------------- " << std::endl;
        std::cout << task->bestNet().format(OutputStyle::NET) << std::endl;
        // std::cout << task->bestNet().format(OutputStyle::TERMINAL) << std::endl;
        std::cout << std::endl;
        std::cout << "Execution Time: " << duration.count() << " milliseconds" << " ---->  merit Value " << task->bestMeritValue() << std::endl;
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
