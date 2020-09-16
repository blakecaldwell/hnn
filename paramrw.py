import numpy as np
from hnn_core import read_params
import fileio as fio

# check if using ongoing inputs
def usingOngoingInputs (params, lty = ['_prox', '_dist']):
  if params is None:
    return False

  try:
    tstop = float(params['tstop'])
  except KeyError:
    return False

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
            print('usingOngoingInputs:',params[dpref[postfix]+k])
            return True
  except: 
    return False
  return False

# return number of evoked inputs (proximal, distal)
# using dictionary d (or if d is a string, first load the dictionary from filename d)
def countEvokedInputs (params):
  nprox = ndist = 0
  if params is not None:
    for k,v in params.items():
      if k.startswith('t_'):
        if k.count('evprox') > 0:
          nprox += 1
        elif k.count('evdist') > 0:
          ndist += 1
  return nprox, ndist

# check if using any evoked inputs 
def usingEvokedInputs (params, lsuffty = ['_evprox_', '_evdist_']):
  nprox,ndist = countEvokedInputs(params)
  if nprox == 0 and ndist == 0:
    return False

  try:
    tstop = float(params['tstop'])
  except KeyError:
    return False

  lsuff = []
  if '_evprox_' in lsuffty:
    for i in range(1,nprox+1,1): lsuff.append('_evprox_'+str(i))
  if '_evdist_' in lsuffty:
    for i in range(1,ndist+1,1): lsuff.append('_evdist_'+str(i))
  for suff in lsuff:
    k = 't' + suff
    if k not in params: continue
    if float(params[k]) > tstop: continue
    k = 'gbar' + suff
    for k1 in params.keys():
      if k1.startswith(k):
        if float(params[k1]) > 0.0: return True
  return False

# check if using any poisson inputs 
def usingPoissonInputs (params):
  if params is None:
    return False

  try:
    tstop = float(params['tstop'])

    if 't0_pois' in params and 'T_pois' in params:
      t0_pois = float(params['t0_pois'])
      if t0_pois > tstop: return False
      T_pois = float(params['T_pois'])
      if t0_pois > T_pois and T_pois != -1.0:
        return False
  except KeyError:
    return False

  for cty in ['L2Pyr', 'L2Basket', 'L5Pyr', 'L5Basket']:
    for sy in ['ampa','nmda']:
      k = cty+'_Pois_A_weight_'+sy
      if k in params:
        if float(params[k]) != 0.0:
          return True

  return False

# check if using any tonic (IClamp) inputs 
def usingTonicInputs (params):
  if params is None:
    return False

  tstop = float(params['tstop'])
  for cty in ['L2Pyr', 'L2Basket', 'L5Pyr', 'L5Basket']:
    k = 'Itonic_A_' + cty + '_soma'
    if k in params:
      amp = float(params[k])
      if amp != 0.0:
        print(k,'amp != 0.0',amp)
        k = 'Itonic_t0_' + cty
        t0,t1 = 0.0,-1.0
        if k in params: t0 = float(params[k])
        k = 'Itonic_T_' + cty
        if k in params: t1 = float(params[k])
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


def consolidate_chunks(inputs):
    # get a list of sorted chunks
    sorted_inputs = sorted(inputs.items(), key=lambda x: x[1]['start'])

    consolidated_chunks = []
    for one_input in sorted_inputs:
        # extract info from sorted list
        input_dict = {'inputs': [one_input[0]],
                      'start': one_input[1]['start'],
                      'end': one_input[1]['end'],
                      'mean': one_input[1]['mean'],
                      'sigma': one_input[1]['sigma'],
                      'opt_start': one_input[1]['opt_start'],
                      'opt_end': one_input[1]['opt_end'],
                      'weights': one_input[1]['weights'],
                      }

        if (len(consolidated_chunks) > 0) and \
            (input_dict['start'] <= consolidated_chunks[-1]['end']):
            # update previous chunk
            consolidated_chunks[-1]['inputs'].extend(input_dict['inputs'])
            consolidated_chunks[-1]['end'] = input_dict['end']
            consolidated_chunks[-1]['opt_end'] = max(consolidated_chunks[-1]['opt_end'], input_dict['opt_end'])
            # average the weights
            consolidated_chunks[-1]['weights'] = (consolidated_chunks[-1]['weights'] + one_input[1]['weights'])/2
        else:
            # new chunk
            consolidated_chunks.append(input_dict)

    return consolidated_chunks

def combine_chunks(input_chunks):
    # Used for creating the opt params of the last step with all inputs

    combined_chunk = {'inputs': [],
                      'opt_start': 0.0,
                      'opt_end': 0.0,
    }

    for evinput in input_chunks:
        combined_chunk['inputs'].extend(evinput['inputs'])
        if evinput['opt_end'] > combined_chunk['opt_end']:
            combined_chunk['opt_end'] = evinput['opt_end']

    # wRMSE with weights of 1's is the same as regular RMSE.
    combined_chunk['weights'] = np.ones(len(input_chunks[-1]['weights']))
    return combined_chunk

def chunk_evinputs(opt_params, sim_tstop, sim_dt):
    import re
    import scipy.stats as stats
    from math import ceil, floor

    num_step = ceil(sim_tstop / sim_dt) + 1
    times = np.linspace(0, sim_tstop, num_step)

    for input_name in opt_params.keys():
        # calculate cdf using start time (minival of optimization range)
        cdf = stats.norm.cdf(times, opt_params[input_name]['start'],
                             opt_params[input_name]['sigma'])
        opt_params[input_name]['cdf'] = cdf.copy()

    for input_name in opt_params.keys():
        opt_params[input_name]['weights'] = opt_params[input_name]['cdf'].copy()

        for other_input in opt_params:
            if input_name == other_input:
                # don't subtract our own cdf(s)
                continue
            if opt_params[other_input]['mean'] < \
               opt_params[input_name]['mean']:
                # check ordering to only use inputs after us
                continue
            else:
                decay_factor = opt_params[input_name]['decay_multiplier']*(opt_params[other_input]['mean'] - \
                                  opt_params[input_name]['mean']) / \
                                  sim_tstop
                opt_params[input_name]['weights'] -= opt_params[other_input]['cdf'] * decay_factor

        # weights should not drop below 0
        opt_params[input_name]['weights'] = np.clip(opt_params[input_name]['weights'], a_min=0, a_max=None)

        # start and stop optimization where the weights are insignificant
        opt_params[input_name]['opt_start'] = min(opt_params[input_name]['start'], times[np.where( opt_params[input_name]['weights'] > 0.01)][0])
        opt_params[input_name]['opt_end'] = max(opt_params[input_name]['end'], times[np.where( opt_params[input_name]['weights'] > 0.01)][-1])

        # convert to multiples of dt
        opt_params[input_name]['opt_start'] = floor(opt_params[input_name]['opt_start']/sim_dt)*sim_dt
        opt_params[input_name]['opt_end'] = ceil(opt_params[input_name]['opt_end']/sim_dt)*sim_dt

    # combined chunks that have overlapping ranges
    # opt_params is a dict, turn into a list
    input_chunks = consolidate_chunks(opt_params)

    # add one last chunk to the end
    if len(input_chunks) > 1:
        input_chunks.append(combine_chunks(input_chunks))

    return input_chunks

def get_inputs (params):
    import re
    input_list = []

    # first pass through all params to get mu and sigma for each
    for k in params.keys():
        input_mu = re.match('^t_ev(prox|dist)_([0-9]+)', k)
        if input_mu:
            id_str = 'ev' + input_mu.group(1) + '_' + input_mu.group(2)
            input_list.append(id_str)

    return input_list

def trans_input (input_var):
    import re

    input_str = input_var
    input_match = re.match('^ev(prox|dist)_([0-9]+)', input_var)
    if input_match:
        if input_match.group(1) == "prox":
            input_str = 'Proximal ' + input_match.group(2)
        if input_match.group(1) == "dist":
            input_str = 'Distal ' + input_match.group(2)

    return input_str