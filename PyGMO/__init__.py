# Copyright (C) 2007, 2008 by Francesco Biscani
# bluescarni@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

from _PyGMO import *

class __vector_init:
	def build(self, iterable, t):
		if iterable == None:
			return
		if not getattr(iterable, '__iter__', False):
			raise TypeError, 'I need an iterable object for initialisation.'
		for i in iterable: self.append(t(i))


class vector_double(_PyGMO.__base_vector_double,__vector_init):
	def __init__(self, iterable = None):
		super(type(self), self).__init__()
		self.build(iterable,float)

class vector_size_t(_PyGMO.__base_vector_size_t,__vector_init):
	def __init__(self, iterable = None):
		super(type(self), self).__init__()
		self.build(iterable,int)
