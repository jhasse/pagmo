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

#ifndef GORADIATOR_H_INCLUDED
#define GORADIATOR_H_INCLUDED

#include <vector>
#include "rad_objfun.h"
#include "GOproblem.h"

//***********************************************************************************
//Nanostructured Radiator Problems....
//***********************************************************************************

class radiatorProb : public GOProblem {
public:
	radiatorProb(int dim);
	virtual double objfun(const std::vector<double>& x) { return radiator(x); }
};	//end class radiatorProb

#endif // GORADIATOR_H_INCLUDED
