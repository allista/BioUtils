'''
Created on Jan 30, 2016

@author: allis
'''

import re, os

replace = r'''(#\n)*# Copyright \(C\) 201\d Allis Tauri <allista@gmail\.com>
# 
# \w* is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# \(at your option\) any later version.
# 
# \w* is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE\.
# See the GNU General Public License for more details\.
# 
# You should have received a copy of the GNU General Public License along
# with this program. +If not, see <http://www\.gnu\.org/licenses/>\.'''

lre = re.compile(replace, re.MULTILINE)

if __name__ == '__main__':
    for dirpath, dirnames, filenames in os.walk(os.path.abspath(os.curdir)):
        for filename in filenames:
            if filename == os.path.basename(__file__): continue
            curfile = os.path.join(dirpath, filename)
            with open(curfile) as inp:
                filetext = inp.read()
                if not lre.search(filetext): continue
                newtext = lre.sub('', filetext)
            with open(curfile, 'w') as out: out.write(newtext)
            print '%s edited' % filename
    print 'Done'
