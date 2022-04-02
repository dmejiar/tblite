# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.

if get_option('api') or get_option('python')
  subdir('api')
endif
subdir('basis')

if get_option('ilp64')
  subdir('blas64')
else
  subdir('blas')
endif

subdir('classical')
subdir('container')
subdir('context')
subdir('coulomb')
subdir('data')
subdir('disp')
subdir('external')
subdir('fit')
subdir('integral')
subdir('io')

if get_option('ilp64')
  subdir('lapack64')
else
  subdir('lapack')
endif

subdir('ncoord')
subdir('output')
subdir('param')
subdir('repulsion')
subdir('scf')
subdir('wavefunction')
subdir('xtb')

srcs += files(
  'adjlist.f90',
  'blas.f90',
  'container.f90',
  'context.f90',
  'cutoff.f90',
  'disp.f90',
  'lapack.f90',
  'ncoord.f90',
  'os.F90',
  'param.f90',
  'repulsion.f90',
  'scf.f90',
  'spin.f90',
  'timer.f90',
  'toml.f90',
  'version.f90',
  'wignerseitz.f90',
)
