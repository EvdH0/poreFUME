#!/bin/bash


#This is a basal integration test. It assumes the user correctly installed the dependencies in INSTALL.md, ran install.sh and env.sh in the current shell. This script will download a set of 75 sequences. All the sequences belong to pacbio barcode 01. A ~40 MB set of raw fast5 files will also be download, since fast5 files are needed for nanopolish. 
#poreFUME is fired up, check whether this matches your local settings, such as cores and directory paths.
#Finally, the output (CARD annotation) is compared to the pre-computed CARD annotation (in test/data) and should be similar.


curdir=`pwd`
echo $curdir

cd test
cd data
wget http://www.student.dtu.dk/~evand/poreFUME_data/testSet75.tar.gz
tar -zxvf testSet75.tar.gz

cd $curdir
cd inputData
wget http://www.student.dtu.dk/~evand/poreFUME_data/testSet75.fasta

cd $curdir
python poreFUME.py inputData/testSet75.fasta inputData/pb_39.fasta --PacBioLegacyBarcode --cores 8 --pathCARD=inputData/n.fasta.protein.homolog.fasta --pathNanocorrect=/home/ubuntu/poreFUME/nanocorrect/ --pathRawreads=/home/ubuntu/poreFUMEtestInstall/poreFUME/test/data/testSet75 --overwriteNanocorrect --pathNanopolish=/home/ubuntu/poreFUMEtestInstall/poreFUME/nanopolish/ --overwriteNanopolish --overwriteDemux --overwriteCARD

DIFF=$(diff output/annotation/testSet75/testSet75.afterNP.annotated.csv test/data/testSet75.afterNP.annotated.csv) 
if [ "$DIFF" != "" ] 
then
    echo "Integration test faile, the output in output/annotation/testSet75/testSet75.afterNP.annotated.csv is not what is expected"
else
    echo "Integration test passed!"
fi
 
