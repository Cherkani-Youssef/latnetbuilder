// This file is part of Nettice Builder.
//
// Copyright (C) 2012-2018  Pierre L'Ecuyer and Universite de Montreal
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef NETBUILDER__PARSER__COMMAND_LINE_H
#define NETBUILDER__PARSER__COMMAND_LINE_H

#include "netbuilder/Types.h"
#include "netbuilder/NetConstructionTraits.h"
#include "netbuilder/Task/Task.h"
#include "netbuilder/FigureOfMerit/FigureOfMerit.h"

#include "latbuilder/SizeParam.h"

namespace NetBuilder { namespace Parser {

/**
 * Collection of arguments required to construct a Search instance.
 */
// template <NetBuilder::NetConstruction , NetBuilder::EmbeddingType>
// struct CommandLine;

/**
 * Specialization of CommandLine for ordinary nets.
 */
template <NetConstruction NC, EmbeddingType ET>
struct CommandLine {
   std::string s_explorationMethod;
   std::string s_designParameter;
   std::string s_size;
   std::string s_dimension;
   std::vector<std::string> s_figures;
   std::string s_figureCombiner;
   std::string s_combiner;
   std::string s_verbose;
   
   bool m_earlyAbort;
   LatBuilder::SizeParam<LatBuilder::LatticeType::DIGITAL, ET> m_sizeParam;
   typename NetConstructionTraits<NC>::DesignParameter m_designParameter;
   Combiner m_combiner;
   Dimension m_dimension;
   std::unique_ptr<FigureOfMerit::FigureOfMerit> m_figure;
   int m_verbose;

   std::unique_ptr<Task::Task> parse();
};

}
}
#include "latbuilder/Parser/SizeParam.h"

#include "netbuilder/Parser/DesignParameterParser.h"
#include "netbuilder/Parser/FigureParser.h"
#include "netbuilder/Parser/ExplorationMethodParser.h"

namespace NetBuilder { namespace Parser {
template <NetConstruction NC, EmbeddingType ET>
std::unique_ptr<NetBuilder::Task::Task>
CommandLine<NC, ET>::parse()
{
      namespace lbp = LatBuilder::Parser;
      m_sizeParam = lbp::SizeParam<LatBuilder::LatticeType::DIGITAL, ET>::parse(s_size);
      m_designParameter = DesignParameterParser<NC,ET>::parse(*this);
      m_dimension = boost::lexical_cast<Dimension>(s_dimension);
      m_verbose = boost::lexical_cast<int>(s_verbose);
      m_figure = FigureParser<NC, ET>::parse(*this); // m_combiner initialized as a side effect
      return ExplorationMethodParser<NC, ET>::parse(*this); // as a side effect, m_figure has been moved to task
}


}}

#endif
