/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

// 21/09/2008: Initial version by Dario Izzo.

#ifndef PAGMO_SOLVERSTHREADS_H
#define PAGMO_SOLVERSTHREADS_H

#include <boost/thread/condition_variable.hpp>
#include <boost/thread/mutex.hpp>

#include "../atomic_counters/atomic_counters.h"
#include "GOproblem.h"
#include "population.h"

//Here we define the parameters needed to instanciate a thread. These contain
//datas that are algorithm specific, but also data that are needed for all aglorithms (LB,UB,objfun,mutex etc.)

struct threadParam{
	//Thread unique ID
	unsigned int threadID;

	//Solvers Data
	int NP;
	int generations;
	//DE
	int strategy;
	double F,CR;
	//PSO
	double omega,eta1,eta2,vcoeff;
	int nswarms;
	//GA
	double M,CRsga;
	int insert_best;
	//SA-AN
	double Ts,Tf;
	//pointers giving access to global resources
	GOProblem* problem;

	atomic_counter_size_t *isActive;
	boost::mutex *TPmutex;
	boost::condition_variable *exit;
	Population *Ptr_pop;
	std::ofstream *Ptr_log;
	uint32_t randomSeed;
};

//Here we define the protoypes for each type of thread we may want to open
void *DEthread(void *data);
void *PSOthread(void *data);
void *MPSOthread(void *data);
void *SGAthread(void *data);
void *ASAthread(void *data);
#endif
