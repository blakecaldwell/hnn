# Human Neocortical Neurosolver (HNN)

## About

The **Human Neocortical Neurosolver (HNN)** is an open-source neural modeling tool designed to help
researchers/clinicians interpret human brain imaging data. HNN presents a convenient GUI to an
anatomically and biophysically detailed model of human thalamocortical brain circuits, which
makes it easier to generate and evaluate hypotheses on the mechanistic origin of signals measured
with MEG/EEG or intracranial ECoG. A unique feature of HNN's model is that it accounts for the
biophysics generating the primary electric currents underlying such data, so simulation results
are directly comparable to source localized data (nano-Ampere-meters); this enables precise
tuning of model parameters to match characteristics of recorded signals.

We are integrating the circuit-level modeling with the minimum-norm-estimate (MNE) source
localization software (https://martinos.org/mne/stable/index.html), so researchers can compute
MEG/EEG source estimates and test hypotheses on
the circuit origin of their data in one software package. Our goal is to design HNN to be useful
to researchers with no formal computational neural modeling or coding experience.

For more information visit [https://hnn.brown.edu](https://hnn.brown.edu) . There, we describe the use of HNN in studying the
circuit-level origin of some of the most commonly measured MEG/EEG and ECoG signal: event related
potentials (ERPs) and low frequency rhythms (alpha/beta/gamma).

## Installation

Please follow the links on our [installation page](installer) to find instructions for your operating system.

## Questions

For questions, comments/feedback, or troubleshooting information please contact
us at hnneurosolver@gmail.com, and review our user forum at
[https://www.neuron.yale.edu/phpBB/viewforum.php?f=46](https://www.neuron.yale.edu/phpBB/viewforum.php?f=46) .

## References

To cite the HNN software please use the following references:
 [eLife 2020;9:e51214 DOI: 10.7554/eLife.51214 ](https://doi.org/10.7554/eLife.51214)
 and
 [![DOI](https://zenodo.org/badge/128077928.svg)](https://zenodo.org/badge/latestdoi/128077928)
 
[![Build Status](https://travis-ci.com/jonescompneurolab/hnn.svg?branch=master)](https://travis-ci.com/jonescompneurolab/hnn)