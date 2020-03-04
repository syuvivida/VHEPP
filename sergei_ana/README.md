## Install the HepSim software toolkit

The HepSim data files are listed in https://atlaswww.hep.anl.gov/hepsim/index.php

### Java 11 installation in Ubuntu
The basic steps are modified from this link:
https://thishosting.rocks/install-java-ubuntu/

First download jdk tarball from (eg. jdk-11.0.6_linux-x64_bin.tar.gz):
https://www.oracle.com/java/technologies/javase-jdk11-downloads.html

Note, you need to have an Oracle account to download the file.

Then, follow the steps below.
```
sudo apt-get update && sudo apt-get upgrade
sudo add-apt-repository ppa:linuxuprising/java
sudo apt-get update
sudo mkdir -p /var/cache/oracle-jdk11-installer-local
sudo cp jdk-11.0.6_linux-x64_bin.tar.gz /var/cache/oracle-jdk11-installer-local/.
sudo apt-get install oracle-java11-installer-local
java -version
```

### Download hs-tool kit
Follow the instruction from this link:
https://atlaswww.hep.anl.gov/hepsim/doc/doku.php?id=hepsim:quick

```
bash
wget http://atlaswww.hep.anl.gov/hepsim/soft/hs-toolkit.tgz -O - | tar -xz;
```
Comment out the following lines in hs-toolkit/setup.sh

```
    if [[ "$version" < "1.7" ]]; then
        echo "Java version is less than 1.7. Please update it."; exit 0;
    fi

```
Now you could do the setup:
```
source hs-toolkit/setup.sh
```
Try a few hs commands to see if your setup is OK.


## HepSim Docker image in Ubuntu
### Install singularity
Follow the basic steps in this link: https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps

#### Install system dependencies
```
sudo apt-get update && sudo apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup
```

#### Install Go
```
bash
export VERSION=1.13 OS=linux ARCH=amd64 
wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz 
sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz 
rm go$VERSION.$OS-$ARCH.tar.gz
echo 'export PATH=/usr/local/go/bin:$PATH' >> ~/.bashrc &&  source ~/.bashrc
```

#### Singularity
```
bash
export VERSION=3.5.2 
wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz 
tar -xzf singularity-${VERSION}.tar.gz
cd singularity

./mconfig && make -C builddir && sudo make -C builddir install
singularity help
```

### Download image file (this could take one hour to half a day) and use singularity
```
wget http://atlaswww.hep.anl.gov/hepsim/soft/centos7hepsim.img
```

Try a few commands in https://atlaswww.hep.anl.gov/hepsim/doc/doku.php?id=hepsim:dev_singularity

```
singularity exec centos7hepsim.img bash -l
```

### Analyzing data

#### Download Sergei's ana example

```
tar xvzf ana_jets_time.tgz
cd ana_jets_time
```
Modify the path of centos7hepsim.img in ana_jets_time/msetup.sh, then

```
bash
source msetup.sh
make
```

#### Add missing steps and missing files
Follow the instruction in ana_jets_time/README up to step 4).

Note, now you are in the directory of ana_jets_time.

Move all the data files:
```
mkdir -p data/pgun_pi
mv data/pgun*slcio data/pgun_pi/.
```
Download detector geometry file for the corresponding input dataset
```
wget https://atlaswww.hep.anl.gov/hepsim/soft/detectors/sifcch7.zip && unzip sifcch7.zip sifcch7.pandora && rm -rf sifcch7.zip
mkdir -p data/rfull009_sifcch7/sifcch7/
mv sifcch7.pandora data/rfull009_sifcch7/sifcch7/.
```

Now you need to modify the data file names in ana_jets_time/A_RUN. Replace the following lines:
```
DAT="pgun_pi10gev"
```

with the energy you want to run on or add extra lines to run on multiple energies in one job
```
DAT="pgun_pi100gev"
ls -1 data/pgun_pi/$DAT* > data.in
./ana
mv root/output.root root/$DAT.root
```


### Run ana example on the downloaded data

Continue from step 5) of ana_jets_time/README and come back to the directory where you store centos7hepsim.img 

```
singularity exec centos7hepsim.img bash -l
cd ana_jets_time
./A_RUN
```

### Study the root files in ana_jets_time/root directory