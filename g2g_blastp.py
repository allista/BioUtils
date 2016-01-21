#!/usr/bin/python
# coding=utf-8
#
# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# degen_primer is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# degen_primer is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
Created on Dec 19, 2015

@author: Allis Tauri <allista@gmail.com>
'''

import os
import signal
from time import sleep

from BioUtils.NCBI import BlastCLI
from reportlab.lib import colors

from BioUtils.Tools.tmpStorage import roDict, clean_tmp_files, shelf_result

_pid = -1
abort_event = None
def sig_handler(signal, frame):
    if _pid != os.getpid(): return
    print('\nAborting. This may take some time '
          'as not all operations could be stopped immediately.\n')
    abort_event.set(); sleep(0.1)
    clean_tmp_files()
#end def

if __name__ == '__main__':
    from multiprocessing import Event
    from BioUtils.Tools.Output import user_message
    from BioUtils.SeqUtils import load_files
    _pid = os.getpid()
    #setup signal handler
    signal.signal(signal.SIGINT,  sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)
    signal.signal(signal.SIGQUIT, sig_handler)
    
    if True:
    #    from DegenPrimer import MultiprocessingBase
    #    MultiprocessingBase.cpu_count = 1
        abort_event = Event()
        lb = BlastCLI(abort_event)
        
        with user_message('Loading genomes...', '\n'):
            genomes_dir = u'/home/allis/Dropbox/Science/Микра/Thermococcus/sequence/GenBank/Thermococcus'
            genome_names = ['Thermococcus_barophilus_Ch5-complete.gb', 
                            'Thermococcus_onnurineus_NA1-complete-genome.gb',
                            'Thermococcus_sp._ES1.gb',
                            'Thermococcus-DS1-preliminary.gb'] 
            genomes = load_files(abort_event, [os.path.join(genomes_dir, f) for f in genome_names], 'gb') 
        
        ref = genomes[0]
        subj = genomes[1:]
        
        @shelf_result
        def g2g2shelf():
            return lb.g2g_blastp(ref, subj, 11, features_of_interest=[{'ugene_name': 'FC-full'}, {'ugene_name': 'COC-full'}])
            
        g2g_res = '/tmp/DP-PCR-N_KAEs'
        if not os.path.isfile(g2g_res):
            g2g_res = g2g2shelf()
            print g2g_res
        if g2g_res:
            with roDict(g2g_res) as db:
                results = db['result'] 
            if results:
                lb.g2g_to_csv('g2g_test.csv', ref, subj, results)
        print 'Done.'
            
    import numpy as np
    import pandas
    import matplotlib.pyplot as plt
    from matplotlib import cm, rc
    
    df = pandas.read_csv('g2g_test.csv')
    df = df.dropna(subset=('Thermococcus_onnurineus_NA1_percent', 'Thermococcus_sp._ES1_percent'))
    ratio = 15/18.0
    rc('text', usetex=True)
    rc('font', size=14)
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    fig, ax = plt.subplots(figsize=(12,12*ratio))
    ax.set_xlabel(r'\textit{Thermococcus paralvinellae} \textsf{ES1}', fontsize=18)
    ax.set_ylabel(r'\textit{Thermococcus onnurineus} \textsf{NA1}', fontsize=18)
    plt.hist2d(np.array(df['Thermococcus_sp._ES1_percent']), 
               np.array(df.Thermococcus_onnurineus_NA1_percent),
               range=[[20,100], [20,100]],
               bins=80, cmap=cm.get_cmap('Blues'))
    plt.colorbar()
    plt.plot([20,100], [20,100], color='lightgrey', linestyle='--')
    plt.plot([50,50], [20,100], color='lightgrey', linestyle='--')
    plt.plot([20,100], [50,50], color='lightgrey', linestyle='--')
    colors = ['darkorange']*2+['darkred']+['red']*6+['orange','darkgreen','darkviolet']+['black']*7
    genes = df.Thermococcus_barophilus_Ch5_gene.isin(range(1737,1756))
    genex = df['Thermococcus_sp._ES1_percent'][genes]
    geney = df.Thermococcus_onnurineus_NA1_percent[genes]
    plt.scatter(x=genex, y=geney, s=60, color=colors[::-1])
    labels = ['fdhA', '4Fe-4S', 
              'MbhH', 'MbhH\'', 'MbhH\'', 'MbhH\'\'', 'MbhM', 'Mbh(K+L)', 'MbhN', 
              'MbhJ', 'MbhX', 'FocA',
              'MbhB', 'MbhC', 'MbhD', 'Mbh(E+F)', 'MbhG', 'MbhA', 'MbhH\'\'\'']
    for label, x, y, c, gid in zip(labels[::-1], genex, geney, colors[::-1], df.Thermococcus_barophilus_Ch5_gene[genes]):
        print gid, label, c
        plt.annotate(
            label, 
            xy = (x, y), xytext = (-3, 3),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            color=c, size=14
#            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
#            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0')
            )
    plt.tight_layout()
    plt.savefig('test.png', dpi=100)
    plt.savefig('test.svg')
    plt.savefig('test.eps', rasterize=True)
    plt.show()