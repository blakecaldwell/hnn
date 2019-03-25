# HNN native install (Mac OS)

This method will run HNN without using virtualization, meaning the GUI may feel more responsive and simulations may run slightly faster. However, the procedure is a set of steps that the user must follow, and there is a possibility that differences in the base environment may require additional troubleshooting. Thus, it is best suited for advanced users. For the recommended Docker-based installation please see the instructions below.
  - Alternative: [Docker install instructions](README.md)

## Prerequisite 1: Xcode Command Line Tools

The Xcode Command Line Tools package includes utilities for compiling code from the terminal window command line (gcc, make, git, etc.). This is needed for compiling mod files in NEURON. To install the package, type the following from a terminal window:
```
xcode-select --install
```


## Prerequisite 2: XQuartz
1. Download the installer image (version 2.7.11 tested): https://www.xquartz.org/
2. Run the XQuartz.pkg installer within the image, granting privileges when requested.
3. Start the XQuartz application. An "X" icon will appear in the taskbar along with a terminal, signaling that XQuartz is waiting for connections. You can minimize the terminal, but do not close it.

From the command line:
```
cd /tmp/
curl https://dl.bintray.com/xquartz/downloads/XQuartz-2.7.11.dmg -o XQuartz-2.7.11.dmg
hdiutil attach /tmp/XQuartz-2.7.11.dmg
sudo installer -verbose -pkg /Volumes/XQuartz-2.7.11/XQuartz.pkg -target /
hdiutil detach /Volumes/XQuartz-2.7.11
rm /tmp/XQuartz-2.7.11.dmg
open /Applications/Utilities/XQuartz.app
```

## Prerequisite 3: Miniconda (Python 3)

1. Run the commands below from a terminal window (as a regular user) or download and install miniconda from the link: https://conda.io/en/latest/miniconda.html. This will create a python environment isolated from other installations on the system (e.g. those installed using homebrew). You could use `brew install python3` if you wish (has been tested with HNN), but this guide will cover the miniconda version.
    ```
    cd /tmp/
    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    sh ./Miniconda3-latest-MacOSX-x86_64.sh -b
    rm /tmp/Miniconda3-latest-MacOSX-x86_64.sh
    ```

## Prerequisite 4: NEURON

1. Install [NEURON](https://neuron.yale.edu/neuron/) from the terminal:
    ```
    cd /tmp/
    curl https://neuron.yale.edu/ftp/neuron/versions/v7.6/nrn-7.6.x86_64-osx.pkg -o nrn-7.6.x86_64-osx.pkg
    sudo installer -verbose -pkg /tmp/nrn-7.6.x86_64-osx.pkg -allowUntrusted -target /
    ```
2. You will be asked about setting PATH variable. Say 'No' to both prompts.

  <img src="install_pngs/neuron_path.png" width="600" />

3. Afterward, you will be presented with a confirmation message that NEURON has been installed. Click 'Continue'

  <img src="install_pngs/neuron_continue.png" width="600" />

## Prepare the Python environment

1. Open a new terminal winodw and run the command below to create a conda environment with the Python prerequisites for HNN.

    ```
    conda create -n hnn python=3.6 mpi4py pyqtgraph pyopengl matplotlib scipy
    ```
2. Run the following from a terminal window (bash shell) to set up the environment variables for future reactivations

    ```
    source activate hnn
    mkdir -p ${CONDA_PREFIX}/etc/conda/activate.d
    mkdir -p ${CONDA_PREFIX}/etc/conda/deactivate.d
    ACTIVATE="${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh"
    DEACTIVATE="${CONDA_PREFIX}/etc/conda/deactivate.d/env_vars.sh"
    echo "export CONDA_OLDPYTHONPATH=\$PYTHONPATH" > $ACTIVATE
    echo "export PYTHONPATH=/Applications/NEURON-7.6/nrn/lib/python:\$PYTHONPATH" >> $ACTIVATE
    echo "export PATH=/Applications/NEURON-7.6/nrn/x86_64/bin:\$PATH" >> $ACTIVATE
    echo "export LD_LIBRARY_PATH=\${CONDA_PREFIX}/lib" >> $ACTIVATE
    echo "unset NRN_PYLIB" >> $ACTIVATE
    echo "unset LD_LIBRARY_PATH" > $DEACTIVATE
    echo "export PYTHONPATH=\$CONDA_OLDPYTHONPATH" >> $DEACTIVATE
    ```

## Reboot your system
Please reboot your system before proceeding. A reboot is really needed after installing Xquartz. The environment variables set above will also be set for all terminal windows after the reboot.

## Compile HNN
1. After rebooting, open a new terminal and clone the repository:
    ```
    git clone https://github.com/jonescompneurolab/hnn.git
    ```
2. Now enter the hnn directory and compile HNN's mod files for NEURON. This is where Xcode Command Line Tools are needed.
    ```
    cd hnn
    make
    ```

## Run the HNN model
1. Start the HNN GUI from a terminal window:
    ```
    source activate hnn
    python hnn.py hnn.cfg
    ```
2. The HNN GUI should appear and you should now be able to run the tutorials at https://hnn.brown.edu/index.php/tutorials/
3. When you run simulations for the first time, the following dialog boxes may pop-up and ask you for permission to allow connections through the firewall. Saying 'Deny' is fine since simulations will just run locally on your Mac.

<img src="install_pngs/nrniv_firewall.png" width="400" />

<img src="install_pngs/orterun_firewall.png" width="400" />

# Troubleshooting

For Mac OS specific issues: please see the [Mac OS troubleshooting page](troubleshooting.md)

If you run into other issues with the installation, please [open an issue on our GitHub](https://github.com/jonescompneurolab/hnn/issues). Our team monitors these issues and will investigate possible fixes.

For other HNN software issues, please visit the [HNN bullentin board](https://www.neuron.yale.edu/phpBB/viewforum.php?f=46)