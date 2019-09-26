import numpy as np
from hnn_core import Params
import fileio as fio

# get dict of ':' separated params from fn; ignore lines starting with #
def quickreadprm (fn):
  return Params(fn)

# get dict of ':' separated params from fn; ignore lines starting with #
def quickgetprm (fn,k,ty):
  d = quickreadprm(fn)
  return ty(d[k])

# check if using ongoing inputs
def usingOngoingInputs (fn, lty = ['_prox', '_dist']):
  params = quickreadprm(fn)
  tstop = float(params['tstop'])
  dpref = {'_prox':'input_prox_A_','_dist':'input_dist_A_'}
  try:
    for postfix in lty:
      if float(params['t0_input'+postfix])<= tstop and \
         float(params['tstop_input'+postfix])>=float(params['t0_input'+postfix]) and \
         float(params['f_input'+postfix])>0.:
        for k in ['weight_L2Pyr_ampa','weight_L2Pyr_nmda',\
                  'weight_L5Pyr_ampa','weight_L5Pyr_nmda',\
                  'weight_inh_ampa','weight_inh_nmda']:
          if float(params[dpref[postfix]+k])>0.:
            print('usingOngoingInputs:',d[dpref[postfix]+k])
            return True
  except: 
    return False
  return False

# return number of evoked inputs (proximal, distal)
# using dictionary d (or if d is a string, first load the dictionary from filename d)
def countEvokedInputs (d):
  if type(d) == str: d = quickreadprm(d)
  nprox = ndist = 0
  for k,v in d.items():
    if k.startswith('t_'):
      if k.count('evprox') > 0:
        nprox += 1
      elif k.count('evdist') > 0:
        ndist += 1
  return nprox, ndist

# check if using any evoked inputs 
def usingEvokedInputs (d, lsuffty = ['_evprox_', '_evdist_']):
  if type(d) == str: d = quickreadprm(d)
  nprox,ndist = countEvokedInputs(d)
  tstop = float(d['tstop']) 
  lsuff = []
  if '_evprox_' in lsuffty:
    for i in range(1,nprox+1,1): lsuff.append('_evprox_'+str(i))
  if '_evdist_' in lsuffty:
    for i in range(1,ndist+1,1): lsuff.append('_evdist_'+str(i))
  for suff in lsuff:
    k = 't' + suff
    if k not in d: continue
    if float(d[k]) > tstop: continue
    k = 'gbar' + suff
    for k1 in d.keys():
      if k1.startswith(k):
        if float(d[k1]) > 0.0: return True
  return False

# check if using any poisson inputs 
def usingPoissonInputs (d):
  if type(d)==str: d = quickreadprm(d)
  tstop = float(d['tstop'])
  if 't0_pois' in d and 'T_pois' in d:
    t0_pois = float(d['t0_pois'])
    if t0_pois > tstop: return False
    T_pois = float(d['T_pois'])
    if t0_pois > T_pois and T_pois != -1.0:
      return False
  for cty in ['L2Pyr', 'L2Basket', 'L5Pyr', 'L5Basket']:
    for sy in ['ampa','nmda']:
      k = cty+'_Pois_A_weight_'+sy
      if k in d:
        if float(d[k]) != 0.0: return True
  return False

# check if using any tonic (IClamp) inputs 
def usingTonicInputs (d):
  if type(d)==str: d = quickreadprm(d)
  tstop = float(d['tstop'])
  for cty in ['L2Pyr', 'L2Basket', 'L5Pyr', 'L5Basket']:
    k = 'Itonic_A_' + cty + '_soma'
    if k in d:
      amp = float(d[k])
      if amp != 0.0:
        print(k,'amp != 0.0',amp)
        k = 'Itonic_t0_' + cty
        t0,t1 = 0.0,-1.0
        if k in d: t0 = float(d[k])
        k = 'Itonic_T_' + cty
        if k in d: t1 = float(d[k])
        if t0 > tstop: continue
        #print('t0:',t0,'t1:',t1)
        if t0 < t1 or t1 == -1.0: return True
  return False

# reads params from a generated txt file and returns gid dict and p dict 
def read_gid_dict (fgid_dict):
    lines = fio.clean_lines(fgid_dict)
    gid_dict = {}
    for line in lines:
        if line.startswith('#'): continue
        keystring, val = line.split(": ")
        key = keystring.strip()
        if val[0] is '[':
            val_range = val[1:-1].split(', ')
            if len(val_range) is 2:
                ind_start = int(val_range[0])
                ind_end = int(val_range[1]) + 1
                gid_dict[key] = np.arange(ind_start, ind_end)
            else:
                gid_dict[key] = np.array([])
    return gid_dict


def write_gid_dict (fgid_dict, gid_dict):
  with open(fgid_dict, 'w') as f:
    pstring = '%26s: '
    # write the gid info first
    for key in gid_dict.keys():
      f.write(pstring % key)
      if len(gid_dict[key]):
        f.write('[%4i, %4i] ' % (gid_dict[key][0], gid_dict[key][-1]))
      else:
        f.write('[]')
      f.write('\n')

# reads params from a generated txt file and returns gid dict and p dict 
def read (fparam):
    lines = fio.clean_lines(fparam)
    p = {}
    gid_dict = {}
    for line in lines:
        if line.startswith('#'): continue
        keystring, val = line.split(": ")
        key = keystring.strip()
        if val[0] is '[':
            val_range = val[1:-1].split(', ')
            if len(val_range) is 2:
                ind_start = int(val_range[0])
                ind_end = int(val_range[1]) + 1
                gid_dict[key] = np.arange(ind_start, ind_end)
            else:
                gid_dict[key] = np.array([])
        else:
            try:
                p[key] = float(val)
            except ValueError:
                p[key] = str(val)
    return gid_dict, p

# write the params to a filename
def write(fparam, p, gid_list):
  """ now sorting
  """
  # sort the items in the dict by key
  # p_sorted = [item for item in p.items()]
  p_keys = [key for key, val in p.items()]
  p_sorted = [(key, p[key]) for key in p_keys]
  # for some reason this is now crashing in python/mpi
  # specifically, lambda sorting in place?
  # p_sorted = [item for item in p.items()]
  # p_sorted.sort(key=lambda x: x[0])
  # open the file for writing
  with open(fparam, 'w') as f:
    pstring = '%26s: '
    # write the gid info first
    for key in gid_list.keys():
      f.write(pstring % key)
      if len(gid_list[key]):
        f.write('[%4i, %4i] ' % (gid_list[key][0], gid_list[key][-1]))
      else:
        f.write('[]')
      f.write('\n')
    # do the params in p_sorted
    for param in p_sorted:
      key, val = param
      f.write(pstring % key)
      if key.startswith('N_'):
        f.write('%i\n' % val)
      else:
        f.write(str(val)+'\n')
