# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# BioUtils is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# indicator_gddccontrol is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Created on Jul 20, 2012

@author: Allis Tauri <allista@gmail.com>
'''

import string
import random


def print_exception(e):
    print '\n%s: error code %s; %s\n' % (type(e).__name__,
                                         hasattr(e, 'errno') and str(e.errno) or 'N/A',  
                                         e.message or 'no message.')
#end def    


def random_text(length):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _unused in xrange(length))
#end def
