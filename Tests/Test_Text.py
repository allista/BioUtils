# coding=utf-8

'''
Created on 2016-01-14

@author: Allis Tauri <allista@gmail.com>
'''

import  cProfile
from BioUtils.Tools.Text import wrap_text, line_by_line

if __name__ == '__main__':
    txt = '''If true, TextWrapper attempts to detect sentence endings and ensure 
    that sentences are always separated by exactly two spaces. This is generally 
    desired for text in a monospaced font. However, the sentence detection    
    algorithm is imperfect: it assumes that a sentence ending consists of a 
    lowercase letter followed by one of '.', '!', or '?', possibly followed by 
    one of '"' or "'", followed by a space. One problem with this is algorithm 
    is that it is unable to detect the difference between “Dr.” in'''
    print wrap_text(txt)
    print '='*80
    
    texts = ['asdfasd  sreydstnsr mywyy    eratg AG RADFG SDFGA lkjoiuguivuasdfhpwoiefjahgaiohghouygmuimjiuh', 
             'qasf; a[r uq[ewjrfasdhfuiah [WERJ AOSF BA;We werqt', 
             'wreyqeuqy a  poiertqprohg aspoei  toeaprt a']
    widths = [20,30,30]
    print line_by_line(texts, widths)
    
    cProfile.run('''for i in xrange(10000): 
    wrap_text(txt)''', 
    'word_wrap.profile')
    print 'Done'