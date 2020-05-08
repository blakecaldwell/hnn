# HNN native install (Windows)

This method will run HNN without using virtualization, meaning the GUI may feel more responsive and simulations may run slightly faster. However, the procedure is a set of steps that the user must follow, and there is a possibility that differences in the base environment may require additional troubleshooting. Thus, it is best suited for advanced users. For the recommended Docker-based installation please see the instructions below.

- Alternative: [Docker install instructions](README.md)

A [PowerShell install script](hnn.ps1) will manage downloading all prerequisites except Microsoft MPI which requires a web browser to download. If the script finds msmpisetup.exe in the Downloads folder, it will take care of installing it.

## Requirements

- A 64-bit OS
- Windows 7 or later. Windows Vista is not supported for lack of multiprocessing support.
- PowerShell version 1.0 or later. If PowerShell is not installed, please follow [this link](https://docs.microsoft.com/en-us/powershell/scripting/install/installing-powershell) for downloading and running the PowerShell installer.

## Run install script

The PowerShell script used below will create a new directory called "hnn" in the place where the command is run from. If you have already cloned a copy of the HNN source code, you can avoid creating this new directory by running the script within the existing source code directory (using the third option below).

1. Run the script from a cmd prompt:

    ```powershell
    @"%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -InputFormat None -ExecutionPolicy Bypass -Command "iex ((New-Object System.Net.WebClient).DownloadString('https://raw.githubusercontent.com/jonescompneurolab/hnn/master/installer/windows/hnn.ps1'))"
    ```

    OR from a powershell prompt:

    ```powershell
    Set-ExecutionPolicy Bypass -Scope Process -Force; iex ((New-Object System.Net.WebClient).DownloadString('https://raw.githubusercontent.com/jonescompneurolab/hnn/master/installer/windows/hnn.ps1'))
    ```

    OR from a local copy (already checked out with git):

    ```powershell
    cd hnn
    powershell.exe -ExecutionPolicy Bypass -File .\installer\windows\hnn.ps1
    ```

   - There will be a permission prompt to install Microsoft MPI and a couple of terminal windows will
open up. There will be a prompt for pressing ENTER after nrnmech.dll has been built
   - If an existing Python 3.X installation isn't found, you should expect that installation will pause for ~5min while installing Miniconda

2. After the script has completed, instructions will be displayed for using the environment either with virtualenv or Miniconda. Open up a new cmd.exe window (not PowerShell) for the environment variables to get set in the session.
3. Run:

    ```powershell
    activate hnn
    cd hnn
    python hnn.py
    ```

4. That will launch the HNN GUI. You should now be able to run the tutorials at https://hnn.brown.edu/index.php/tutorials/

## Troubleshooting

### Running hnn fails with "Permission denied" for python3

When trying to run simulations in HNN, you might see messages similar to below:

```powershell
Starting simulation (2 cores). . .
Simulation exited with return code 4294967293. Stderr from console:
NEURON -- VERSION 7.6.5 master (f3dad62b) 2019-01-11
Duke, Yale, and the BlueBrain Project -- Copyright 1984-2018
See http://neuron.yale.edu/neuron/credits

C:/nrn/bin/nrnpyenv.sh: line 141: /cygdrive/c/Users/[USERNAME]/AppData/Local/Microsoft/WindowsApps/python3: Permission denied
Python not available
```

This issue occurs after a particular Windows update is applied that inserts a non-functional alias for Python into the PATH environment variable before the functional version we installed. Since these are non-functional aliases, it is fine to remove them. HNN will then be able to locate the correct Python executable.

```powershell
rm c:\Users\[USERNAME]\AppData\Local\Microsoft\WindowsApps\python.exe
rm c:\Users\[USERNAME]\AppData\Local\Microsoft\WindowsApps\python3.exe
```

### Other issues

If you run into other issues with the installation, please [open an issue on our GitHub](https://github.com/jonescompneurolab/hnn/issues). Our team monitors these issues and will investigate possible fixes.

For other HNN software issues, please visit the [HNN bulletin board](https://www.neuron.yale.edu/phpBB/viewforum.php?f=46)
