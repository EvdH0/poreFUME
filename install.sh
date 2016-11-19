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
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.4.0/ncbi-blast-2.4.0+-x64-linux.tar.gz
tar -xzf ncbi-blast-2.4.0+-x64-linux.tar.gz

#Install samtools
cd $curdir
git clone --recursive https://github.com/samtools/htslib.git
cd htslib; git checkout 6bed35a3eaefa3baa2c7e0166ceba442212f166b;make

#probably you need: apt-get install libncurses5-dev
cd $curdir
git clone --recursive https://github.com/samtools/samtools.git
cd samtools; git checkout 897c0027a04501e3ea33d94b5cdeb633d010da8d; make

#Instal bwa
cd $curdir
git clone https://github.com/lh3/bwa.git
cd bwa; git checkout 5961611c358e480110793bbf241523a3cfac049b; make

#Install nanopolish, we used a forked version that has a 'build in' supporting fraction filter
cd $curdir
git clone --recursive https://github.com/EvdH0/nanopolish.git
cd nanopolish; git checkout 04fd9aecbb4ab266350476b957f4abb8ed994d8d; make 

#Install GNU parallel, we don't actually use this for nanopolish, but kept it in as a possiblity
cd $curdir
wget http://ftp.gnu.org/gnu/parallel/parallel-20160922.tar.bz2
tar xvjf parallel-20160922.tar.bz2
cd parallel-20160922/
./configure && make && sudo make install
