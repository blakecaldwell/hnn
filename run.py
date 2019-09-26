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
import newparamrw d
import specfn as specfn
#import pickle
import datetime
from conf import readconf
from L5_pyramidal import L5Pyr
from L2_pyramidal import L2Pyr
from L2_basket import L2Basket
from L5_basket import L5Basket

import os.path as op

###############################################################################
# Let us import hnn_core

import hnn_core
from hnn_core import simulate_dipole, average_dipoles, get_rank, shutdown
from hnn_core import Params, Network

dconf = readconf()
doutf = {}

# TODO: should be moved to Network class in hnn core
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

def savedat (params, dpl, net, spikedata):
  global doutf

  newparamrw.write(doutf['file_param'], params, net.gid_dict)

  dpl.write(doutf['file_dpl_norm'])

  # TODO: this should be moved to Network class within hnn-core
  # write the somatic current to a file
  # for now does not write the total but just L2 somatic and L5 somatic
  X = np.r_[[dpl.t, net.current['L2Pyr_soma'].x, net.current['L5Pyr_soma'].x]].T
  np.savetxt(doutf['file_current'], X, fmt=['%3.3f', '%5.4f', '%5.4f'],
              delimiter='\t')

  # write spikes file
  np.savetxt(doutf['file_spikes'], spikedata, fmt=['%3.2f', '%d'],
              delimiter='\t')

  # if p['save_vsoma']: save_vsoma()

  # for i,elec in enumerate(lelec):
  #   elec.lfpout(fn=doutf['file_lfp'].split('.txt')[0]+'_'+str(i)+'.txt',tvec = t_vec)


def runanalysis (params, dpl, fspec):
  spec_opts = {'type': 'dpl_laminar',
                'f_max': params['f_max_spec'],
                'save_data': 0,
                'runtype': 'parallel',
              }
  specfn.analysis_simp(spec_opts, params, dpl, fspec) # run the spectral analysis

def setupsimdir (params, f_psim):
  simdir = os.path.join(dconf['datdir'], params['sim_prefix'])
  try:
    os.mkdir(simdir)
  except FileExistsError:
    pass

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
  params = Params(f_psim)

  ddir = setupsimdir(params, f_psim) # one directory for all experiments
  # create rotating data files
  setoutfiles(ddir)

  net = Network(params)
  arrangelayers(net) # arrange cells in layers - for visualization purposes

  dpls = simulate_dipole(net, params['N_trials'])
  spikedata = net.allreduce_spikes()

  if get_rank() == 0:
    if dpls == None:
      print("ERR: Failed to start simulate_dipole")
      exit(2)

    avg_dpl = average_dipoles(dpls)

  if params['save_spec_data'] or newparamrw.usingOngoingInputs(f_psim):
    spec = doutf['file_spec']
    runanalysis(params, avg_dpl, spec) # run spectral analysis

  if get_rank() == 0:
    savedat(params, avg_dpl, net, spikedata)

    # below are not updated for hnn-core yet
    # for elec in lelec: print('end; t_vec.size()',t_vec.size(),'elec.lfp_t.size()',elec.lfp_t.size())

    # if params['save_figs']:
    #   savefigs(params) # save output figures

if __name__ == "__main__":
  f_psim = None

  # LFP is not working with hnn-core
  # testLFP = dconf['testlfp']
  # testlaminarLFP = dconf['testlaminarlfp']
  # lelec = [] # list of LFP electrodes

  # reads the specified param file
  for i in range(len(sys.argv)):
    if sys.argv[i].endswith('.json'):
      f_psim = sys.argv[i]
      break
    if sys.argv[i].endswith('.param'):
      f_psim = sys.argv[i]
      break

  if f_psim is None:
    f_psim = os.path.join('param','default.param')

  if os.path.exists(f_psim):
    runsim(f_psim)
  else:
    print("ERR: could not find param file: %s" % os.path.normpath(f_psim))
    exit(1)

  # terminate
  shutdown()