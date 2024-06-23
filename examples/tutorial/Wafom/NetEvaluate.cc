

#include <iostream>
#include <memory>

#include "netbuilder/Types.h"
#include "netbuilder/DigitalNet.h"
#include "netbuilder/FigureOfMerit/WeightedFigureOfMerit.h"

#include "netbuilder/FigureOfMerit/Wafom/ProjMeritFASTWAFOM.h"
#include "netbuilder/FigureOfMerit/Wafom/FastWAFOMProjMerit.h"
#include "netbuilder/FigureOfMerit/Wafom/Wafom.h"
#include "netbuilder/FigureOfMerit/Wafom/FastWafom.h"

#include "Path.h"

#include "netbuilder/Task/Eval.h"

#include "latbuilder/Util.h"

#include "latbuilder/Functor/LookUpTable.h"

#include "netbuilder/Types.h"
#include "netbuilder/Helpers/JoeKuo.h"
#include "netbuilder/FigureOfMerit/TValue.h"
#include "netbuilder/FigureOfMerit/CoordUniformFigureOfMerit.h"

using namespace NetBuilder;
using JoeKuo::createJoeKuoSobolNet;
using namespace NetBuilder::FigureOfMerit;
using namespace NetBuilder::Task;
// namespace plt = matplotlibcpp;

int main(int argc, char **argv)
{
    SET_PATH_TO_LATNETBUILDER_FOR_EXAMPLES();
    uInteger k, w, q;
    Dimension dim;
    double h, factor;

    if (argc != 7)
    {
        // Parse command line arguments
        w = 15;
        k = w;
        q = 3;
        dim = 10;

        h = 1.0;
        factor = 1.0;
    }
    else
    {
        // Parse command line arguments
        w = std::atoi(argv[1]);
        k = std::atoi(argv[2]);
        q = std::atoi(argv[3]);
        dim = std::atoi(argv[4]);
        h = std::atof(argv[5]);
        factor = std::atof(argv[6]);
    }
    int l = w / q;

    {
        LookUpTable table_c(w, q, l, h, factor);

        auto figure = std::make_unique<FastWafom>(q, w, std::move(table_c));

        std::cout << "Computing the Wafom..." << std::endl;

        std::cout << figure->format() << std::endl;
        auto evaluator = figure->evaluator();
        auto net = createJoeKuoSobolNet(dim, k);
        std::cout << "Merit value: " << (*evaluator)(net) << std::endl
                  << std::endl;
    }
}