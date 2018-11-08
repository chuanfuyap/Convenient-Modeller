Bootstrap: docker
From: ubuntu:16.04

%files
dist.zip /opt

%post
apt -y update
apt-get -y install autotools-dev
apt-get -y install autoconf
apt-get -y install default-jre
apt-get -y install unzip
apt-get -y install wget
apt-get -y install build-essential
apt-get -y install libxml2-dev
apt-get -y install libtool
apt-get -y install pkg-config

cd
wget -O sundials-2.4.0.tar.gz https://computation.llnl.gov/projects/sundials/download/sundials-2.4.0.tar.gz
tar -xzf sundials-2.4.0.tar.gz
cd sundials-2.4.0
./configure --enable-shared CFLAGS=-fPIC
make
make install

cd
wget -O libSBML-5.13.0-core-src.tar.gz https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/libSBML-5.13.0-core-src.tar.gz/download
tar -xzf libSBML-5.13.0-core-src.tar.gz
cd libsbml-5.13.0
./configure
make
make install

cd
wget -O SBML_odeSolver.zip https://github.com/raim/SBML_odeSolver/archive/master.zip
unzip SBML_odeSolver.zip
cd SBML_odeSolver-master
autoreconf -i
./configure
make
make install

cd /opt
unzip dist.zip

cd /opt
wget -O soslibJNIC.zip https://github.com/chuanfuyap/soslibJNIC/archive/master.zip
unzip soslibJNIC.zip
cd soslibJNIC-master/
make
cp ./dist/sosLibLinkv3.so /opt/lib

%environment
export LD_LIBRARY_PATH=/usr/local/lib

%runscript
exec java -jar /opt/GRaPe2.jar gui
