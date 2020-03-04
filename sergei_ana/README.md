## Install the HepSim software toolkit

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