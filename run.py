#!/usr/bin/env python
# run.py - primary run function for s1 project
#
# v 1.10.0-py35
# rev 2016-05-01 (SL: removed izip, fixed an nhost bug)
# last major: (SL: toward python3)
# other branch for hnn

import os
import sys
import time
import shutil
import numpy as np
# Cells are defined in other files
import network
import fileio as fio
from newparamrw import usingOngoingInputs
import plotfn as plotfn
import specfn as specfn
import pickle
import datetime
from dipolefn import Dipole
from conf import readconf
from L5_pyramidal import L5Pyr
from L2_pyramidal import L2Pyr
from L2_basket import L2Basket
from L5_basket import L5Basket
from lfp import LFPElectrode

import os.path as op

###############################################################################
# Let us import hnn_core

import hnn_core
from hnn_core import simulate_dipole, Params, Network, shutdown

dconf = readconf()
doutf = {}

# NEEDS to be moved inside hnn-core
# spike write function
# def spikes_write (net, filename_spikes):
#   f = open(filename_spikes,'w')
#   f.close() # first make sure writes to an empty file
#   for rank in range(int(pc.nhost())):
#     # guarantees node order and no competition
#     pc.barrier()
#     if rank == int(pc.id()):
#       # net.spiketimes and net.spikegids are type h.Vector()
#       L = int(net.spikegids.size())
#       with open(filename_spikes, 'a') as file_spikes:
#         for i in range(L):
#           file_spikes.write('%3.2f\t%d\n' % (net.spiketimes.x[i], net.spikegids.x[i]))
#   # let all nodes iterate through loop in which only one rank writes
#   pc.barrier()

# # save somatic voltage of all cells to pkl object
# def save_vsoma ():
#   for host in range(int(pc.nhost())):
#     if host == get_rank():
#       dsoma = net.get_vsoma()
#       messageid = pc.pack(dsoma) # create a message ID and store this value
#       pc.post(host,messageid) # post the message
#   if get_rank()==0:
#     dsomaout = {}
#     for host in range(int(pc.nhost())):
#       pc.take(host)
#       dsoma_node = pc.upkpyobj()
#       for k,v in dsoma_node.items(): dsomaout[k] = v
#     dsomaout['vtime'] = t_vec.to_python()
#     # print('dsomaout.keys():',dsomaout.keys(),'file:',doutf['file_vsoma'])
#     pickle.dump(dsomaout,open(doutf['file_vsoma'],'wb'))

# copies param file into root dsim directory
def copy_paramfile (dsim, f_psim, str_date):
  fout = os.path.join(dsim,f_psim.split(os.path.sep)[-1])
  shutil.copyfile(f_psim,fout)
  # open the new param file and append the date to it
  with open(fout, 'a') as f_param: f_param.write('\nRun_Date: %s' % str_date)

def savedat (p, dpl, net):
  global doutf
  """
  Writing the dipole output is now part of hnn-core.simulate_dipole, but the
  remaining code for writing somatic current, output spikes, and lfp need
  to be updated
  """

  pass
  # write the somatic current to the file
  # for now does not write the total but just L2 somatic and L5 somatic
  # with open(doutf['file_current'], 'w') as fc:
  #   for t, i_L2, i_L5 in zip(t_vec.x, net.current['L2Pyr_soma'].x, net.current['L5Pyr_soma'].x):
  #     fc.write("%03.3f\t" % t)
  #     # fc.write("%5.4f\t" % (i_L2 + i_L5))
  #     fc.write("%5.4f\t" % i_L2)
  #     fc.write("%5.4f\n" % i_L5)
  # write output spikes
  # file_spikes_tmp = fio.file_spike_tmp(dproj)
  # spikes_write(net, file_spikes_tmp)
  # # move the spike file to the spike dir
  # if rank == 0: shutil.move(file_spikes_tmp, doutf['file_spikes'])
  # if p['save_vsoma']: save_vsoma()
  # for i,elec in enumerate(lelec):
  #   elec.lfpout(fn=doutf['file_lfp'].split('.txt')[0]+'_'+str(i)+'.txt',tvec = t_vec)


def runanalysis (prm, fparam, fdpl, fspec):
  """ Needs to be updated for hnn-core """
  #   if get_rank()==0: print("Running spectral analysis...",)
  #   spec_opts = {'type': 'dpl_laminar',
  #                'f_max': prm['f_max_spec'],
  #                'save_data': 0,
  #                'runtype': 'parallel',
  #              }
  #   t_start_analysis = time.time()
  #   specfn.analysis_simp(spec_opts, fparam, fdpl, fspec) # run the spectral analysis
  #   if get_rank()==0 and debug: print("time: %4.4f s" % (time.time() - t_start_analysis))


def setupsimdir (params, f_psim):
  simdir = os.path.join(dconf['datdir'], params['sim_prefix'])
  try:
    os.mkdir(simdir)
  except FileExistsError:
    pass
  str_date = datetime.datetime.now().strftime("%Y-%m-%d")
  #copy_paramfile(simdir, f_psim, str_date)
  return simdir

def getfname (ddir,key,trial=0,ntrial=1):
  datatypes = {'rawspk': ('spk','.txt'),
               'rawdpl': ('rawdpl','.txt'),
               'normdpl': ('dpl','.txt'), # same output name - do not need both raw and normalized dipole - unless debugging
               'rawcurrent': ('i','.txt'),
               'rawspec': ('rawspec','.npz'),
               'rawspeccurrent': ('speci','.npz'),
               'avgdpl': ('dplavg','.txt'),
               'avgspec': ('specavg','.npz'),
               'figavgdpl': ('dplavg','.png'),
               'figavgspec': ('specavg','.png'),
               'figdpl': ('dpl','.png'),
               'figspec': ('spec','.png'),
               'figspk': ('spk','.png'),
               'param': ('param','.txt'),
               'vsoma': ('vsoma','.pkl'),
               'lfp': ('lfp', '.txt')
             }
  if ntrial == 1 or key == 'param': # param file currently identical for all trials
    return os.path.join(ddir,datatypes[key][0]+datatypes[key][1])
  else:
    return os.path.join(ddir,datatypes[key][0] + '_' + str(trial) + datatypes[key][1])
    

# create file names
def setoutfiles (ddir,trial=0,ntrial=1):
  # if get_rank()==0: print('setoutfiles:',trial,ntrial)
  global doutf
  doutf['file_dpl'] = getfname(ddir,'rawdpl',trial,ntrial)
  doutf['file_current'] = getfname(ddir,'rawcurrent',trial,ntrial)
  doutf['file_param'] = getfname(ddir, 'param',trial,ntrial)
  doutf['file_spikes'] = getfname(ddir, 'rawspk',trial,ntrial)
  doutf['file_spec'] = getfname(ddir, 'rawspec',trial,ntrial)
  doutf['filename_debug'] = 'debug.dat'
  doutf['file_dpl_norm'] = getfname(ddir,'normdpl',trial,ntrial)
  doutf['file_vsoma'] = getfname(ddir,'vsoma',trial,ntrial)
  doutf['file_lfp'] = getfname(ddir,'lfp',trial,ntrial)
  #if get_rank()==0: print('doutf:',doutf)
  return doutf

def expandbbox (boxA, boxB):
  return [(min(boxA[i][0],boxB[i][0]),max(boxA[i][1],boxB[i][1]))  for i in range(3)]

def arrangelayers (net):
  # offsets for L2, L5 cells so that L5 below L2 in display
  dyoff = {L2Pyr: 1000, 'L2_pyramidal' : 1000,
           L5Pyr: -1000-149.39990234375, 'L5_pyramidal' : -1000-149.39990234375,
           L2Basket: 1000, 'L2_basket' : 1000,
           L5Basket: -1000-149.39990234375, 'L5_basket' : -1000-149.39990234375}
  for cell in net.cells: cell.translate3d(0,dyoff[cell.celltype],0)
  dbbox = {x:[[1e9,-1e9],[1e9,-1e9],[1e9,-1e9]] for x in dyoff.keys()}
  for cell in net.cells:
    dbbox[cell.celltype] = expandbbox(dbbox[cell.celltype], cell.getbbox())

# All units for time: ms
def runsim (f_psim):
  t0 = time.time() # clock start time

  params = Params(f_psim, json_fmt=False)

  ddir = setupsimdir(params, f_psim) # one directory for all experiments
  # create rotating data files
  setoutfiles(ddir)

  net = Network(params)
  arrangelayers(net) # arrange cells in layers - for visualization purposes

  dpl = simulate_dipole(net)

  shutdown()

  if debug: print("Simulation run time: %4.4f s" % (time.time()-t0))
  if debug: print("Simulation directory is: %s" % ddir.dsim)
  if params['save_spec_data'] or usingOngoingInputs(doutf['file_param']):
    runanalysis(params, doutf['file_param'], doutf['file_dpl_norm'], doutf['file_spec']) # run spectral analysis

  # below is not updated for hnn-core yet

  # savedat(params, dpl, net)
  # for elec in lelec: print('end; t_vec.size()',t_vec.size(),'elec.lfp_t.size()',elec.lfp_t.size())

  # if params['save_figs']:
  #   savefigs(params) # save output figures

if __name__ == "__main__":

  if dconf['dorun']:

    hnn_core_root = op.join(op.dirname(hnn_core.__file__), '..')

    # data directory - ./data
    dproj = dconf['datdir'] # fio.return_data_dir(dconf['datdir'])
    debug = dconf['debug']
    f_psim = ''
    ntrial = 1
    simlength = 0.0

    # LFP is not working with hnn-core
    # testLFP = dconf['testlfp']
    # testlaminarLFP = dconf['testlaminarlfp']
    # lelec = [] # list of LFP electrodes

    # reads the specified param file
    foundprm = False
    for i in range(len(sys.argv)):
      if sys.argv[i].endswith('.param'):
        f_psim = sys.argv[i]
        foundprm = True
        if debug: print('using ',f_psim,' param file.')
      elif sys.argv[i] == 'ntrial' and i+1<len(sys.argv):
        ntrial = int(sys.argv[i+1])
        if ntrial < 1: ntrial = 1
        if debug: print('ntrial:',ntrial)
      elif sys.argv[i] == 'simlength' and i+1<len(sys.argv):
        simlength = float(sys.argv[i+1])
        if debug: print('simlength:',simlength)

    if not foundprm:
      f_psim = os.path.join('param','default.param')
      if debug: print(f_psim)

    runsim(f_psim)
