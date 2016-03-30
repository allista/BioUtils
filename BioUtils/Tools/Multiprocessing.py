# coding=utf-8

'''
Created on Mar 6, 2013

@author: Allis Tauri <allista@gmail.com>
'''

import errno
import os, traceback, signal
import multiprocessing as mp
from multiprocessing.queues import Queue
from time import sleep
from collections import Sequence
from threading import Thread
from Queue import Empty
from copy import deepcopy
from functools import partial

from .UMP import UProcess
from .tmpStorage import clean_tmp_files
from .AbortableBase import AbortableBase, aborted
from .Debug import raise_tb, raise_tb_on_error, estr

cpu_count = mp.cpu_count()

#decorators
def worker(func):
    func = raise_tb_on_error(func)
    def worker(queue, abort_event, *args, **kwargs):
        try:
            result = func(abort_event, *args, **kwargs)
            if aborted(abort_event): queue.cancel_join_thread()
            else: queue.put(result)
        except Exception, e:
            print 'Exception in %s:\n%s' % (str(func), estr(e))
            abort_event.set()
            queue.cancel_join_thread()
        finally:
            queue.put(None)
            queue.close()
    return worker
#end def
    
def worker_method(func):
    def worker_method(self, queue, abort_event, *args, **kwargs):
        worker(partial(func, self))(queue, abort_event, *args, **kwargs)
    return worker_method
#end def
    
def data_mapper(func):
    func = raise_tb_on_error(func)
    def mapper(queue, abort_event, in_queue, data, *args, **kwargs):
        if kwargs.get('copy_data', False): data = deepcopy(data)
        init = kwargs.get('init_args')
        if init: args = init(*args)
        init = kwargs.get('init_data')
        if init: data = init(data)
        while True:
            try: 
                item = in_queue.get(True, 0.1)
                if item is None: break 
                result = func(data[item], *args)
                if aborted(abort_event): 
                    queue.cancel_join_thread() 
                    break
                queue.put((item, result))
            except Empty: continue
            except Exception, e:
                print 'Exception in %s\n%s:' % (str(func), estr(e))
                abort_event.set()
                queue.cancel_join_thread()
                break
        queue.put(None)
        queue.close()
    return mapper
#end def
    
def data_mapper_method(func):
    def mapper_method(self, queue, abort_event, in_queue, data, *args, **kwargs):
        data_mapper(partial(func, self))(queue, abort_event, in_queue, data, *args, **kwargs)
    return mapper_method
#end def
    
def results_assembler(func):
    def assembler(result, *args):
        func(result[0], result[1], *args)
    return assembler
#end def
    
def results_assembler_methd(func):
    def assembler_method(self, result, *args):
        func(self, result[0], result[1], *args)
    return assembler_method
#end def


def even_chunks(work, num_jobs):
    work_len = len(work)
    if work_len <= num_jobs: 
        #if work size is less then number of requested jobs, run less jobs
        chunks = [1 for _n in xrange(work_len)]
    else:
        #distribute work evenly between jobs
        ch  = work_len/num_jobs
        rem = work_len%num_jobs
        if rem: 
            chunks   = [ch+1 for _r in xrange(rem)]
            chunks  += [ch for _n in xrange(num_jobs-rem)]
        else: chunks = [ch for _n in xrange(num_jobs)]
    return chunks
#end def

@results_assembler
def ordered_results_assembler(index, result, output):
    output[index] = result
    
@results_assembler
def unordered_results_assembler(index, result, output):
    output.append(result)

class Job(object):
    __slots__ = ('process', 'queue')
    def __init__(self, p, q):
        self.process = p
        self.queue = q
        
    def get(self, wait, timeout):
        return self.queue.get(wait, timeout)
    
    def start(self): 
        if not self.process.is_alive(): 
            self.process.start()
    
    def is_alive(self): return self.process.is_alive()
    def join(self): self.process.join()
    
    
class Work(Sequence, Thread, AbortableBase):
    def __init__(self, abort_event, jobs=None, timeout=1, daemonic=0, **kwargs):
        Thread.__init__(self)
        AbortableBase.__init__(self, abort_event)
        self.daemon     = True
        self._jobs      = jobs or []
        self._daemonic  = daemonic
        self._timeout   = timeout
        self.assembler  = None
        self._assembler_args = None
        self._counter   = kwargs.get('counter')
        self._launched  = False
    #end defs
    
    def __len__(self): return len(self._jobs)
    
    def __getitem__(self, index): return self._jobs[index]
    
    def add_job(self, job, queue): self._jobs.append(Job(job, queue))

    def set_jobs(self, jobs): self._jobs = jobs

    def start_jobs(self):
        if self._launched: return
        for j in self._jobs: j.start()
        self._launched = True
    #end def

    def start_work(self, worker, work, num_jobs, *args, **kwargs):
        '''work should be and indexable sequence''' 
        wlen = len(work)
        if not wlen: return
        if self._counter is not None: self._counter.set_work(wlen)
        #determine number of jobs to start
        if not num_jobs: num_jobs = cpu_count
        #prepare jobs
        in_queue = Queue(wlen+num_jobs)
        self._jobs = [None]*num_jobs
        for j in xrange(num_jobs):
            queue = Queue()
            job   = UProcess(target=worker, args=(queue, self._abort_event, 
                                                  in_queue, work)+args, kwargs=kwargs)
            job.daemon = self._daemonic
            self._jobs[j] = Job(job,queue)
        self.start_jobs()
        for i in xrange(wlen): in_queue.put(i, False)
        for j in xrange(num_jobs): in_queue.put(None, False)
    #end def
    
    def set_assembler(self, assembler, *args):
        assert hasattr(assembler, '__call__'), 'Assembler should be a callable'
        self.assembler = assembler
        self._assembler_args = args
    
    def assemble(self, assembler, *args):
        self.set_assembler(assembler, *args)
        self.get_nowait()
    
    def get_result(self):
        assert self.assembler is not None, 'Assembler should be set before calling get_results'
        if not self._launched: return
        while self._jobs:
            finished_job = None
            for i,job in enumerate(self._jobs):
                try: 
                    out = job.get(True, self._timeout)
                    if out is not None:
                        self.assembler(out, *self._assembler_args)
                        if self._counter is not None: 
                            self._counter.count()
                    else: 
                        job.join()
                        finished_job = i
                    break
                except Empty:
                    if not job.is_alive():
                        job.join()
                        finished_job = i
                    continue
                except IOError, e:
                    if e.errno == errno.EINTR: continue
                    else:
                        print 'Unhandled IOError in Work.get_results:\n%s' % estr(e)
                        raise
                except Exception, e:
                    print 'Unhandled Exception in Work.get_results:\n%s' % estr(e)
                    raise
            if finished_job is not None:
                del self._jobs[finished_job]
        if self._counter is not None and not self.aborted(): 
            self._counter.done()
    #end def
    get_nowait = Thread.start
    
    def run(self):
        try: self.get_result()
        #Ctrl-C
        #EOF means that target activity was terminated in the Manager process
        except (KeyboardInterrupt, EOFError): return
        #IO code=4 means the same
        except IOError, e:
            if e.errno == errno.EINTR: return
            print 'Unhandled IOError in Work thread:\n%s' % estr(e)
            self._abort_event.set()
        except Exception, e:
            print 'Unhandled Exception in Work thread:\n%s' % estr(e)
            self._abort_event.set()
    #end def
    
    def wait(self):
        if self.is_alive(): self.join()
        return not self.aborted()
#end class        
        

class MultiprocessingBase(AbortableBase):
    '''Base class for classes that use job parallelization through multiprocessing'''

    #aliases
    _Process = staticmethod(UProcess)
    _Queue   = staticmethod(Queue)
    
    worker                  = staticmethod(worker)
    worker_method           = staticmethod(worker_method)
    data_mapper             = staticmethod(data_mapper)
    data_mapper_method      = staticmethod(data_mapper_method)
    results_assembler       = staticmethod(results_assembler)
    results_assembler_methd = staticmethod(results_assembler_methd)
    
    raise_tb                = staticmethod(raise_tb)
    raise_tb_on_error       = staticmethod(raise_tb_on_error)
    
    ordered_results_assembler   = staticmethod(ordered_results_assembler)
    unordered_results_assembler = staticmethod(unordered_results_assembler)
    
    #class body
    def __init__(self, abort_event, daemonic=True):
        AbortableBase.__init__(self, abort_event)
        self._daemonic = daemonic
    #end def


    def Work(self, jobs=None, timeout=1, **kwargs): 
        return Work(self._abort_event, jobs, timeout, self._daemonic, **kwargs)
    
    def start_jobs(self, *work):
        for w in work: w.start_jobs()

    def get_nowait(self, *work):
        for w in work: w.get_nowait()
        
    def wait(self, *work):
        ret = True
        for w in work: ret &= w.wait()
        return ret
    
    def parallelize2(self, timeout, worker, assembler, work, *args, **kwargs):
        w = self.Work(timeout=timeout, **kwargs)
        w.start_work(worker, work, None, *args, **kwargs)
        w.assemble(assembler)
        if not w.wait(): return False
        return True
    
    def parallelize(self, timeout, worker, work, *args, **kwargs):
        #prepare and start jobs
        w = self.Work(timeout=timeout, **kwargs)
        w.start_work(worker, work, None, *args, **kwargs)
        #allocate result container, get result
        result = [None]*len(work)
        w.assemble(ordered_results_assembler, result)
        #if aborted, return None
        if not w.wait(): return None
        return result
    #end def

    def parallelize_work(self, timeout, func, data, *args, **kwargs):
        if not data: return None
        if len(data) == 1: return [func(data[0], *args)]
        worker = MultiprocessingBase.data_mapper(func)
        return self.parallelize(timeout, worker, data, *args, **kwargs)
    #end def    
    
    def parallelize_functions(self, timeout, funcs, *args, **kwargs):
        if not funcs: return None
        if len(funcs) == 1: return [funcs[0](*args)]
        @MultiprocessingBase.data_mapper
        def worker(func, *args): return func(*args)
        return self.parallelize(timeout, worker, funcs, *args, **kwargs)
    #end def
    
    def parallelize_both(self, timeout, funcs, data, *args, **kwargs):
        if not data or not funcs: return None
        assert len(funcs) == len(data), 'Number of functions must be equal to the data length'
        if len(data) == 1: return [funcs[0](data[0], *args)]
        @MultiprocessingBase.data_mapper
        def worker(func_and_data, *args):
            func, item = func_and_data
            return func(item, *args)
        return self.parallelize(timeout, worker, zip(funcs,data), *args, **kwargs)
    #end def
#end class


#shortcuts
def parallelize_work(abort_event, daemonic, timeout, func, data, *args, **kwargs):
    '''Parallel map implementation with MultiprocessingBase'''
    return MultiprocessingBase(abort_event, daemonic).parallelize_work(timeout, func, data, *args, **kwargs)

def parallelize_functions(abort_event, daemonic, timeout, funcs, *args, **kwargs):
    '''Execute functions from funcs list in parallel using MultiprocessingBase.'''
    return MultiprocessingBase(abort_event, daemonic).parallelize_functions(timeout, funcs, *args, **kwargs)

def parallelize_both(abort_event, daemonic, timeout, funcs, data, *args, **kwargs):
    '''Execute func_i(data_i, *args) for func_i, data_i in zip(funcs, data) in parallel using MultiprocessingBase.'''
    return MultiprocessingBase(abort_event, daemonic).parallelize_both(timeout, funcs, data, *args, **kwargs)


#decorators and managers
from .tmpStorage import shelf_result
from .UMP import FuncManager
from multiprocessing import Manager
import sys
import argparse

Parallelizer = FuncManager('Parallelizer', 
                           (parallelize_work,
                            parallelize_functions,
                            parallelize_both,
                            shelf_result(parallelize_work),
                            shelf_result(parallelize_functions),
                            shelf_result(parallelize_both)))

class MPMain(object):
    description = 'Description is not filled in a subclass'
    
    def __init__(self, pid=os.getpid(), run=False):
        self.pid = pid
        self.mgr = Manager()
        self.abort_event = self.mgr.Event()
        self.argparser = None
        self.args = None
        signal.signal(signal.SIGINT,  self.sig_handler)
        signal.signal(signal.SIGTERM, self.sig_handler)
        signal.signal(signal.SIGQUIT, self.sig_handler)
        if run: self()
        
    def __del__(self):
        try: self.mgr.stop()
        except: pass
        
    def sig_handler(self, signal, frame):
        if self.pid != os.getpid(): return
        print('\nAborting. This may take some time '
              'as not all operations could be stopped immediately.\n')
        self.abort_event.set(); sleep(0.1)
        clean_tmp_files()
    
    def _main(self): pass
    
    def argument(self, *args, **kwargs):
        if self.argparser is None:
            self.argparser = argparse.ArgumentParser(self.description)
        self.argparser.add_argument(*args, **kwargs)
        
    def parse_args(self):
        self.args = self.argparser.parse_args()
    
    def __call__(self, sys_exit=True, *args, **kwargs):
        try: ret = self._main()
        except SystemExit, e:
            if sys_exit: sys.exit(e.code)
            else: return e.code
        except Exception, e:
            self.abort_event.set()
            print 'Unhandled Exception in MPMain:\n%s' % estr(e) 
            traceback.print_exc()
            if sys_exit: sys.exit(1)
            else: return 1
        if sys_exit: sys.exit(ret or 0)
        else: return 0
#end class