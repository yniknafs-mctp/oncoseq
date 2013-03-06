'''
AssemblyLine: transcriptome meta-assembly from RNA-Seq

Copyright (C) 2012,2013 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

def main():
    setup(name="oncoseq",
          version='0.3.1',
          description="high-throughput sequence data analysis in cancer",
          long_description=__doc__,
          author = "Matthew Iyer",
          author_email = "mkiyer@umich.edu",
          license="GPL3",
          platforms="Linux",
          url="http://oncoseq.googlecode.com",
          packages=['oncoseq',
                    'oncoseq.rnaseq',
                    'oncoseq.rnaseq.lib',
                    'oncoseq.rnaseq.pipeline',
                    'oncoseq.rnaseq.utils'],
          cmdclass= {'build_ext': build_ext})

if __name__ == '__main__': 
    main()