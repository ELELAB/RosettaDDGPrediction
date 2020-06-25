#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
#    Script containing unit tests for run_ddg_protocols.py (to be run
#    with pytest).
#
#    Copyright (C) 2020 Valentina Sora <sora.valentina1@gmail.com>
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program. 
#    If not, see <http://www.gnu.org/licenses/>.

import os
import os.path
import pytest

# append ../src to PYTHONPATH
import sys, os
sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/../src"))

import run_ddg_protocols as rddg

class TestRunDdgProtocols:

    _testfilespath = os.path.join(os.getcwd(), "test_files")
    _pdbfile = "stability.pdb"
    _reslistfile = "reslistfile.txt"
    _listfile = "listfile_stability.txt"
    _ncaalistfile = "ncaalistfile.txt"
    _scorefile = "scorefile.sc"
    _reslist = ["A", "C", "D", "E"]
    _mutlist = [([("A", "L", "52", "F"), ("A", "F", "51", "X[DALA]")],), \
                ([("A", "F", "118", "W")],)]

    @pytest.fixture()
    def path(self):
        return os.getcwd()

    @pytest.fixture()
    def pdbfile(self):
        return os.path.join(self._testfilespath, self._pdbfile)

    @pytest.fixture()
    def reslistfile(self):
        return os.path.join(self._testfilespath, self._reslistfile)

    @pytest.fixture()
    def listfile(self):
        return os.path.join(self._testfilespath, self._listfile)

    @pytest.fixture()
    def ncaalistfile(self):
        return os.path.join(self._testfilespath, self._ncaalistfile)

    @pytest.fixture()
    def scorefile(self):
        return os.path.join(self._testfilespath, self._scorefile)

    @pytest.fixture()
    def reslist(self):
        return self._reslist

    @pytest.fixture()
    def mutlist(self):
        return self._mutlist


    def test_get_abspath(self, path):
        return rddg.get_abspath(path)

    def test_get_reslist(self, reslistfile):
        return rddg.get_reslist(reslistfile)

    def test_get_mutlist(self, listfile, reslist):
        mutlist1 = rddg.get_mutlist(listfile, None, None)
        mutlist2 = rddg.get_mutlist(listfile, reslist, None)
        mutlist3 = rddg.get_mutlist(listfile, None, 2)
        mutlist4 = rddg.get_mutlist(listfile, reslist, 2)
        return [mutlist1, mutlist2, mutlist3, mutlist4]

    def test_get_pose_mutlist(self, pdbfile, mutlist):
        return rddg.get_pose_mutlist(pdbfile, mutlist)

    def test_get_ncaa(self, ncaalistfile):
        return rddg.get_ncaa(ncaalistfile)

    def test_get_struct_lowscore(self, scorefile):
        return rddg.get_struct_lowscore(scorefile)