# based on https://github.com/kenneth59715/nsg-rest-client/blob/master/nsg.nopassword.ipynb
# This works with python 3, with requests module installed
# use port 8443 for production, 8444 for test
# register at https://www.nsgportal.org/reg/reg.php for username and password

import os
import requests
import xml.etree.ElementTree
import time
import sys
import re
import zipfile
import tarfile
import glob
from conf import dconf

debug = dconf['debug']

def getuserpass ():
  f = open('nsgr.txt')
  l = f.readlines()
  CRA_USER = l[0].strip()
  PASSWORD = l[1].strip() # #'changeme'
  f.close()
  return CRA_USER,PASSWORD

CRA_USER,PASSWORD = getuserpass() # this will be collected from the HNN GUI later on

# for production version:
# log in at https://nsgr.sdsc.edu:8443/restusers/login.action
# Tool names can be found at Developer->Documentation (Tools: How to Configure Specific Tools)
# create a new application at Developer->Application Management (Create New Application)
# save the Application Key for use in REST requests

KEY = 'HNN-418776D750A84FC28A19D5EF1C7B4933'
TOOL = 'SINGULARITY_HNN_TG' 
URL = 'https://nsgr.sdsc.edu:8443/cipresrest/v1' # for production version

def createpayload (paramf, ntrial, tstop):
  # returns dictionary of parameters for the NSG job
  payload = {'metadata.statusEmail' : 'true'} 
  payload['vparam.runtime_'] = 0.1 # 0.5
  payload['vparam.filename_'] = 'run.py'
  payload['vparam.cmdlineopts_'] = '-nohomeout -paramf ' + os.path.join('param',paramf) + ' ' + str(ntrial)
  payload['vparam.number_nodes_'] = 1
  payload['tool'] = TOOL
  return payload

#
def prepinputzip (fout='test.zip'):
  """ prepares input zip file for NSGR; file contains all py,mod,param,cfg
  files needed to run the simulation """
  try:
    if debug: print('Preparing NSGR input zip file...',zippath)
    fp = zipfile.ZipFile(fout, "w")
    lglob = ['*.py','mod/*.mod','*.cfg','param/*.param','res/*.png','Makefile']
    for glb in lglob:
      for name in glob.glob(glb):
        if debug: print('adding:',os.path.realpath(name))      
        if name.endswith('.mod'):
          fp.write(name, 'hnn/mod/'+os.path.basename(name), zipfile.ZIP_DEFLATED)
        elif name.endswith('.param'):
          fp.write(name,'hnn/param/'+os.path.basename(name),zipfile.ZIP_DEFLATED)
        elif name.endswith('.png'):
          fp.write(name,'hnn/res/'+os.path.basename(name),zipfile.ZIP_DEFLATED)        
        else:
          fp.write(name, 'hnn/'+os.path.basename(name), zipfile.ZIP_DEFLATED)
    fp.close()
    return True
  except:
    print('prepinputzip ERR: could not prepare input zip file',fout,'for NSGR.')
    return False
    
#
def untar (fname):
  # extract contents of tar gz file to current directory
  tar = tarfile.open(fname)
  tar.extractall()
  tar.close()
  print("Extracted",fname," in Current Directory.")

#
def procoutputtar (fname='output.tar.gz'):
  """ process HNN NSGR output tar file, saving simulation data
  and param file to appropriate directories """
  try:
    tar = tarfile.open(fname)
    for member in tar.getmembers():
      if member.isreg():  # skip if not a file (e.g. directory)
        f = member.name
        if f.count('data')>0:
          lp = f.split(os.path.sep)
          member.name = os.path.basename(member.name) # remove the path by resetting it
          tar.extract(member,os.path.join('data',lp[-2])) # extract to data subdir
          if f.endswith('.param'):
            tar.extract(member,'param') # extract to param subdir
    tar.close()
    if debug: print("Extracted",fname)
    return True
  except:
    print('procoutputtar ERR: Could not extract contents of ',fname)
    return False

def runjobNSGR (paramf='default.param', ntrial=1, tstop=710.0):
  """ run a simulation job on NSG using Restful interface; first prepares input zip
  file, then submits job and waits for it to finish, finally downloads simulation output
  data and extracts it to appropriate location """

  try:

    payload = createpayload(paramf,ntrial,tstop)
    headers = {'cipres-appkey' : KEY} # application KEY
    zippath = os.path.realpath('inputfile.zip')

    if not prepinputzip(zippath):
      print('runjobNSGR ERR: could not prepare NSGR input zip file',zippath)
      return False

    files = {'input.infile_' : open(zippath,'rb')} # input zip file with code to run

    r = requests.post('{}/job/{}'.format(URL, CRA_USER), auth=(CRA_USER, PASSWORD), data=payload, headers=headers, files=files)
    #print(r.text)
    root = xml.etree.ElementTree.fromstring(r.text)

    # sys.stderr.write("%s\n" % r.text)
    sys.stderr.write("%s\n" % r.url)

    for child in root:
      if child.tag == 'resultsUri':
        for urlchild in child:
          if urlchild.tag == 'url':
            outputuri = urlchild.text
      if child.tag == 'selfUri':
        for urlchild in child:
          if urlchild.tag == 'url':
            selfuri = urlchild.text

    #print(outputuri,file=sys.stderr)
    sys.stderr.write("%s\n" % outputuri)
    #print(selfuri,file=sys.stderr)
    sys.stderr.write("%s\n" % selfuri)

    #print('Waiting for job to complete',file=sys.stderr)
    sys.stderr.write('Waiting for job to complete\n')
    jobdone = False
    while not jobdone:
      r = requests.get(selfuri, auth=(CRA_USER, PASSWORD), headers=headers)
      #print(r.text)
      root = xml.etree.ElementTree.fromstring(r.text)
      for child in root:
        if child.tag == 'terminalStage':
          jobstatus = child.text
          if jobstatus == 'false':
            time.sleep(5)
            #print('.',file=sys.stderr,end='')
            sys.stderr.write('.')
          else:
            jobdone = True
            #print('',file=sys.stderr,end='\n')
            sys.stderr.write('\n')
            break

    #print('Job completion detected, getting download URIs...',file=sys.stderr)
    sys.stderr.write('Job completion detected, getting download URIs...')

    r = requests.get(outputuri,
                     headers= headers, auth=(CRA_USER, PASSWORD))
    #print(r.text)
    globaldownloadurilist = []
    root = xml.etree.ElementTree.fromstring(r.text)
    for child in root:
      if child.tag == 'jobfiles':
        for jobchild in child:
          if jobchild.tag == 'jobfile':
            for downloadchild in jobchild:
              if downloadchild.tag == 'downloadUri':
                for attchild in downloadchild:
                  if attchild.tag == 'url':
                    #print(attchild.text)
                    sys.stdout.write(attchild.text)
                    globaldownloadurilist.append(attchild.text)

    #print('Download complete.  Run the next cell.',file=sys.stderr)
    sys.stderr.write('NSG download complete.\n')

    #submitoutput.show()
    #print(submitoutput.stdout)
    #print(globaldownloadurilist)
    globaloutputdict = {}
    for downloaduri in globaldownloadurilist:
      r = requests.get(downloaduri, auth=(CRA_USER, PASSWORD), headers=headers)
      #print(r.text)
      globaloutputdict[downloaduri] = r.text

    #http://stackoverflow.com/questions/31804799/how-to-get-pdf-filename-with-python-requests
    for downloaduri in globaldownloadurilist:
      r = requests.get(downloaduri, auth=(CRA_USER, PASSWORD), headers=headers)
      sys.stderr.write("%s\n" % r.headers)
      d = r.headers['content-disposition']
      fname_list = re.findall("filename=(.+)", d)
      for fname in fname_list:
        if debug: sys.stderr.write("%s\n" % fname)

    # download all output files
    for downloaduri in globaldownloadurilist:
      r = requests.get(downloaduri, auth=(CRA_USER, PASSWORD), headers=headers)
      #sys.stderr.write("%s\n" % r.json())
      #r.content
      d = r.headers['content-disposition']
      filename_list = re.findall('filename=(.+)', d)
      for filename in filename_list:
        #http://docs.python-requests.org/en/master/user/quickstart/#raw-response-content
        with open(filename, 'wb') as fd:
          for chunk in r.iter_content():
            fd.write(chunk)

    # get a list of jobs for user and app key, and terminalStage status
    r = requests.get("%s/job/%s" % (URL,CRA_USER), auth=(CRA_USER, PASSWORD), headers=headers)
    #print(r.text)

    ldeluri = [] # list of jobs to delete
    root = xml.etree.ElementTree.fromstring(r.text)
    for child in root:
      if child.tag == 'jobs':
        for jobchild in child:
          if jobchild.tag == 'jobstatus':
            for statuschild in jobchild:
              if statuschild.tag == 'selfUri':
                for selfchild in statuschild:
                  if selfchild.tag == 'url':
                    #print(child)
                    joburi = selfchild.text
                    jobr = requests.get(joburi, auth=(CRA_USER, PASSWORD), headers=headers)
                    jobroot = xml.etree.ElementTree.fromstring(jobr.text)
                    for jobrchild in jobroot:
                      if jobrchild.tag == 'terminalStage':
                        jobstatus = jobrchild.text
                        sys.stdout.write("job url: %s status terminalStage: %s\n" % (joburi,jobstatus))
                        ldeluri.append(joburi)

    # get information for a single job, print out raw XML, need to set joburi according to above list
    # delete an old job, need to set joburi
    for joburi in ldeluri:
      if debug: print('deleting old job with joburi = ',joburi)
      #joburi = 'https://nsgr.sdsc.edu:8443/cipresrest/v1/job/kenneth/NGBW-JOB-NEURON73_TG-220F7B3C7EE84BC3ADD87346E933ED5E'
      r = requests.get(joburi, headers= headers, auth=(CRA_USER, PASSWORD))
      #print(r.text)
      r = requests.delete(joburi, auth=(CRA_USER, PASSWORD), headers=headers)
      if debug: sys.stderr.write("%s\n" % r.text)

    if not procoutputtar('output.tar.gz'):
      print('runjobNSGR ERR: could not extract simulation output data.')
      return False

    return True

  except:
    print('runjobNSGR unhandled exception!')
    return False
