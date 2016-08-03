#!/bin/bash
curdir=`pwd`
echo $curdir

#For the python dependencies pleace check INSTALL.md 


#Nanocorrect install
git clone https://github.com/jts/nanocorrect.git
ln -s nanocorrect/poa-blosum80.mat
cd nanocorrect; git checkout 0cc9da028156a14892ed592163f647822fe21792;

#Most ugly file update, we make a second copy of nanocorrect that omits the coverage of 3 requirement
cd $curdir
sed '88s/.*/    min_coverage = 1/' nanocorrect/nanocorrect.py > nanocorrect/nanocorrectC0.py

#POA install
cd $curdir
wget http://downloads.sourceforge.net/project/poamsa/poamsa/2.0/poaV2.tar.gz
tar -xzf poaV2.tar.gz
cd poaV2; make CFLAGS='-O3 -g -DUSE_WEIGHTED_LINKS -DUSE_PROJECT_HEADER -I.' poa
cd $curdir
ln -s poaV2/poa

#DALIGNER install
cd $curdir
git clone https://github.com/thegenemyers/DALIGNER.git
cd DALIGNER; git checkout 549da77b91395dd; make

#DAZZ_DB install
cd $curdir
git clone https://github.com/thegenemyers/DAZZ_DB
cd DAZZ_DB; git checkout 8cb2f29c4011a2c2; make


#install blast
cd $curdir
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.4.0+-x64-linux.tar.gz
tar -xzf ncbi-blast-2.4.0+-x64-linux.tar.gz
