// This file is part of Lattice Builder.
//
// Copyright (C) 2012-2016  Pierre L'Ecuyer and Universite de Montreal
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

#ifndef NETBUILDER__TASK__BASE_TASK_H
#define NETBUILDER__TASK__BASE_TASK_H

#include <ostream>
#include <memory>

#include "netbuilder/DigitalNet.h"

namespace NetBuilder { namespace Task {

/**
 * Base base class for all tasks.
 */
class BaseTask {
public:
   virtual ~BaseTask() {}

   /**
    * Executes the task.
    */
    virtual void execute() = 0;

    virtual const DigitalNet& netOutput() const = 0;

    virtual Real meritValueOutput() const = 0;

protected:
   //virtual void format(std::ostream& os) const = 0;
   //friend std::ostream& operator<<(std::ostream&, const Task&);
};

/**
 * Formats and outputs \c task to \c os.
 */
//inline std::ostream& operator<<(std::ostream& os, const Task& task)
//{ task.format(os); return os; }
//

}}

#endif
