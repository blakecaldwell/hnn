# Installing HNN on CentOS (Docker install)

This guide describes installing HNN on CentOS using Docker. This method will automatically download the HNN Docker container image when HNN is started for the first time. If you would prefer to install HNN without Docker, please see the instructions below.
  - Alternative: [Native install instructions (advanced users)](native_install.md)

**Note: please verify that your OS supports Docker CE**
   - CentOS 7 and later versions are supported
   - See [OS requirements](https://docs.docker.com/install/linux/docker-ce/centos/#os-requirements)
   - It may be possible to install older versions of docker meant for earlier CentOS versions. We recommend following the currently supported procedure from Docker, but you may find a version of Docker in the [EPEL repositories](https://fedoraproject.org/wiki/EPEL) for CentOS 6, which would would work for your OS.

## Prerequisite: install Docker
* Open a bash terminal and follow [Docker's CentOS install instructions](https://docs.docker.com/install/linux/docker-ce/centos/) to install the docker-ce RPM packages from Docker's official repository
* Add your user to the docker group to avoid having to run docker commands with "sudo"
    ```
    sudo usermod -a -G docker [username]
    ```
* Log out and back in for the group change to take effect

## Prerequisite: install Docker Compose
* Open a bash terminal and run these commands (from [Docker Compose installation](https://docs.docker.com/compose/install/)):
    ```
    sudo curl -L "https://github.com/docker/compose/releases/download/1.23.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
    sudo chmod +x /usr/local/bin/docker-compose
    docker-compose --version
    ```

## Start HNN

1. Check that Docker is running properly by typing the following in a new terminal window.
    ```
    docker info
    ```
2. Clone or download the [HNN repo](https://github.com/jonescompneurolab/hnn). If you already have a previous version of the repository, bring it up to date with the command `git pull origin master` instead of the `git clone` command below.
    ```
    git clone https://github.com/jonescompneurolab/hnn.git
    cd hnn/installer/docker
    ```
4. For some Linux systems it may be necessary to explicitly authorize the docker container to use the X server display for the GUI. The command below may not always be necessary, but we recommend it to be safe.
    ```
    xhost +local:docker
    ```
5. Start the Docker container. Note: the jonescompneurolab/hnn docker image will be downloaded from Docker Hub (about 1.5 GB). Docker-compose starts a docker container based on the specification file docker-compose.yml and "up" starts the containers in that file and "-d" starts the docker containers in the background.
    ```
    docker-compose up -d
    ```    
6. The HNN GUI should show up and you should now be able to run the tutorials at https://hnn.brown.edu/index.php/tutorials/
   * A directory called "hnn" exists both inside the container (at /home/hnn_user/hnn) and outside (in the directory where step 3 was run) that can be used to share files between the container and your host OS.
   * If you run into problems starting the Docker container or the GUI is not displaying, please see the [Docker troubleshooting section](../docker/README.md#Troubleshooting)
   * If you closed the HNN GUI, and would like to restart it, run the following:
      ```
      docker-compose restart
      ```
7. **NOTE:** You may want run commands or edit files within the container. To access a command prompt in the container, use [`docker exec`](https://docs.docker.com/engine/reference/commandline/exec/):
    ```
    C:\Users\myuser>docker exec -ti docker_hnn_1 bash
    hnn_user@054ba0c64625:/home/hnn_user$
    ```

    If you'd like to be able to copy files from the host OS without using the shared directory, you do so directly with [`docker cp`](https://docs.docker.com/engine/reference/commandline/cp/).

## Uninstalling HNN

1. If you want to just remove the container and 1.5 GB HNN image, run these commands from a terminal window:
    ```
    docker rm -f docker_hnn_1
    docker rmi jonescompneurolab/hnn
    ```
2. To continue and remove Docker, follow these instructions: [Uninstall Docker CE](https://docs.docker.com/install/linux/docker-ce/centos/#uninstall-docker-ce)


# Troubleshooting

If you run into other issues with the installation, please [open an issue on our GitHub](https://github.com/jonescompneurolab/hnn/issues). Our team monitors these issues and will be able to suggest possible fixes.

For other HNN software issues, please visit the [HNN bullentin board](https://www.neuron.yale.edu/phpBB/viewforum.php?f=46)