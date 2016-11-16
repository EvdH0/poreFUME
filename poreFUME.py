from Bio import SeqIO
import pandas as pd
import logging
import sys
import os
import Bio as bio
import operator
import numpy as np
import pandas as pd
import pickle
import math
import errno
from types import *

from time import sleep

from random import randint

#from shutil import rmtree

import shutil

from shutil import copyfile

from Barcode import Barcode
import argparse

import multiprocessing
from subprocess import Popen, PIPE




logger = logging.getLogger()



def main():
    """
    Main function that calls subroutines.
    Layout:
    Setup command line parser
    Setup logging
    
    1a. Call demux()
    1b. Call demuxCollect()
    2.  Call nanocorrect()
    3.  Call nanopolish()
    4.  Call annotateCARD()
    done
    
    All the steps can be turned on and off using --skipXXX this allows flexible analysis and intermidate exits.
    
    1a Demux
    uses the Smith-Waterman implementation which is very slow on this dataset. Relevant parameters are --barcodeEdge , --barcodeThreshold and --match, --mismatch, --gapopen, --gapextend. The returned sequences will be reverse complemented such that the read always starts with the forward primer. The collected reads are stored in mysample.AFTERBC.fasta
    2. nanocorrect uses nanocorrect by Loman,Simpson,Quick and is also slow. Relevant parameters are --pathNanocorrect. The collected reads are stored as mysample.AFTERNC1.fasta (first round of nanocorrect) and mysample.AFTERNC2.fasta (second round of nanocorrect without the coverage requirement)
    3. nanopolish uses nanopolish by Loman,Simpson,Quick. 
    4. annotateCARD find high scoring segmens in the result list, which is not extremly fast. --annotateAll will invoke the annotation of the start and intermidiate files (ie. the afterBC and afterNC1 )
    
    
    
    
    This  is an example of initial file structure
    
    inputData/
        -yourBarcodeData.fasta
        -yourNanoporeData.fasta
    
    
    poreFUME will create
    
    
    output/
        barcode/
            yoursample/
                {n}.yoursample.afterBC.fasta
        etc.
    
    Final results the user is interested in are in output/annotation/xxx where
    
    """
    
    ###
    ### Setup the command line parser
    ###
    parser = argparse.ArgumentParser()
    parser.add_argument("fileONTreads", help="path to FASTA where the (2D) nanopore reads are stored",
                        type=str)
    parser.add_argument("fileBarcodes", help="path to FASTA where the barcodes are stored, format should be ie F_34 for forward and R_34 for reverse barcode",
                        type=str)


    parser.add_argument("--PacBioLegacyBarcode",help="the pacbio_barcodes_paired.fasta file has first digist as 4 instead of 04, turning this option on will fix this",
                        action="store_true")
    
    parser.add_argument("--verbose",help="switch the logging from INFO to DEBUG",
                        action="store_true")
    
    parser.add_argument("--overwriteDemux",help="overwrite results in the output/barcode/runid directory if they exist",
                        action="store_true")


    
    parser.add_argument("--overwriteNanocorrect",help="overwrite the results in the output/nanocorrect/runid directory if the exist",
                                                                action="store_true")

    parser.add_argument("--overwriteNanopolish",help="overwrite the results in the output/nanopolish/runid directory, if the exist",
                                                                action="store_true")
                    

    
    parser.add_argument("--overwriteCARD",help="overwrite the results in the output/annotation/runid directory if the exist",
                                                                action="store_true")
    
    parser.add_argument("--skipDemux",help="Skip the barcode demux step and proceed with nanocorrect, cannot be used with overwrite. Assumes the output/barcode/ and output/ directory are populated accordingly",
                                            action="store_true")
    parser.add_argument("--skipDemuxCollect",help="will skip the demux it self and go to collection based on the pickle",
                                                                action="store_true")
    
    parser.add_argument("--skipNanocorrect",help="Skip the nanocorrect step.",
                                            action="store_true")

    parser.add_argument("--skipNanopolish",help="Skip the nanocorrect step.",
                                            action="store_true")

    parser.add_argument("--skipCARD",help="Skip the CARD annotation",
                                            action="store_true")
                    
    
    parser.add_argument("--match",help="Score for match in alignment (default: %(default)s)",nargs='?', default=2.7,type=float)
    parser.add_argument("--mismatch",help="Score for mis-match in alignment (default: %(default)s)",nargs='?', default=-4.5,type=float)
    parser.add_argument("--gapopen",help="Score for gap-open in alignment (default: %(default)s)",nargs='?', default=-4.7,type=float)
    parser.add_argument("--gapextend",help="Score for gap-extend in alignment (default: %(default)s)",nargs='?', default=-1.6,type=float)
    parser.add_argument("--cores",help="Amount of args.cores to use for multiprocessing (default: %(default)s)",nargs='?', default=1,type=int)
    parser.add_argument("--barcodeThreshold",help="Minimum score for a barcode pair to pass (default: %(default)s)",nargs='?', default=58,type=int) #58 was used on the 'lib A set', 54 on porecamp?
    parser.add_argument("--barcodeEdge",help="Maximum amount of bp from the edge of a read to look for a barcode. (default: %(default)s)",nargs='?', default=60,type=int) #60 was used on the 'lib A set' , 120 on the lib B since it had a different experimental ligation protocol
    

    parser.add_argument("--pathNanocorrect",help="Set the path to the nanocorrect files (default: %(default)s)",nargs='?', default="/Users/evand/Downloads/testnanocorrect/nanocorrect/",type=str) #make sure it has an extra copy of nanocorrect with a coverage of 0 in it
    parser.add_argument("--pathNanopolish",help="Set the path to the nanopolish files (default: %(default)s)",nargs='?', default="/Users/evand/Downloads/nanopolish/nanopolish/",type=str) #location of nanopolish
    parser.add_argument("--pathBWA",help="Set the path to BWA (default: %(default)s)",nargs='?', default="/Users/evand/Downloads/nanopolish/bwa",type=str) #location of nanopolish

    parser.add_argument("--pathRawreads",help="Set the path to the raw reads (.fast5 files), nanopolish needs this. As a hint, this should be the absolute path to which the last part of the header on the poretools produced fasta file referes to. poreFUME will make a symlink to the directory containing the .fast5 files.",nargs='?', default="inputData/NB6",type=str) #location of fast5 files as refered to in the datafiles created by poreTools. See nanopolish docs for more info

    parser.add_argument("--pathCARD",help="Set the path to CARD fasta file (default: %(default)s)",nargs='?', default="inputData/n.fasta.protein.homolog.fasta",type=str) 
    
    parser.add_argument("--annotateAll",help="By default only the  final (demuxed and two times corrected) dataset is annotated, however by turning on this option all the files, raw, after demux, after 1st round of correction, after 2nd round of correction are annotated. This obviously takes longer.",
                                            action="store_true")
    parser.add_argument("--minCoverage",help="sequences will only be nanopolish'ed if they have a coverage that is higher than this threshold. (default: %(default)s)",nargs='?', default=30,type=int) #Jared suggested using at least 30x coverage



    
    args = parser.parse_args()



    #baseFileName is used throughout poreFUME to refer to the specific runname
    try:
        baseFileName =  os.path.splitext(os.path.basename(args.fileONTreads))[0]
    except:
        ValueError('Invalid file name (fileONTreads)')

    
    if not os.path.exists('output'): #it it does not exist
        os.makedirs('output') #create it

    
    ####################
    ### Setup logging ##
    ####################
    
    # create file handler which logs even debug messages
    fh = logging.FileHandler('output/' + baseFileName + '.info.log')
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s -  %(levelname)s - %(message)s')
    fh.setFormatter(formatter)

    
    fh2 = logging.FileHandler('output/' + baseFileName + '.debug.log')
    fh2.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s -  %(levelname)s - %(message)s')
    fh2.setFormatter(formatter)


    
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    if args.verbose: #Switch the root log level
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    
    logging.basicConfig(stream=sys.stdout, level=loglevel,format='%(asctime)s - %(name)s -  %(levelname)s - %(message)s',handlers=[fh,ch])
    
    logging.getLogger('').addHandler(fh)
    logging.getLogger('').addHandler(fh2)







    
    ##############
    #Test inputs #
    ##############
    

    #check if blast is in the path
    if not cmdExists('makeblastdb'): 
        raise RuntimeError('makeblastdb is not avialable in PATH')
    else:
        logger.info('makeblastdb found in path')
    
    if not cmdExists('blastn'):
        raise RuntimeError('blastn is not avialable in PATH')
    else:
        logger.info('blastn found in path')
    
    
    if not cmdExists('poa'):
        raise RuntimeError('poa is not avialable in PATH, install http://sourceforge.net/projects/poamsa/')
    else:
        logger.info('poa found in path')
        
        
    if not cmdExists('LAcat'):
        raise RuntimeError('LAcat is not avialable in PATH, install https://github.com/thegenemyers/DALIGNER and https://github.com/thegenemyers/DAZZ_DB')
    else:
        logger.info('LAcat found in path')
        
    

    
    if not os.path.isfile(args.fileONTreads):
        raise IOError('fileONTreads does not exist',args.fileONTreads)

    
    if not os.path.isfile(args.fileBarcodes):
        raise IOError('fileBarcodes does not exist',args.fileBarcodes)
    
    #Barcode dir
    if not os.path.exists(getBarcodeDir(baseFileName)): #it it does not exist
        os.makedirs(getBarcodeDir(baseFileName)) #create it

    

    
    
    
    if not args.skipNanocorrect:
        if not os.path.exists(args.pathNanocorrect):
            raise IOError('the pathNanocorrect is not valid!')

    ###start putting correct paths inplace
    
    if not args.skipNanocorrect: #Do a check for clean directory now instead of when demux is done
        ##Direcotry handeling
        nanocorrectDir = getNanocorrectDir(baseFileName)
        if not os.path.exists(nanocorrectDir): #it it does not exist
            os.makedirs(nanocorrectDir) #create it
            logger.info(nanocorrectDir + ' did not exist, created it')
        else: #it exists
            logger.info(nanocorrectDir + ' already exists')
            if args.overwriteNanocorrect: #if we are sure to overwrite the existing barcode output
                shutil.rmtree(nanocorrectDir)
                os.makedirs(nanocorrectDir) #create it
                logger.info(nanocorrectDir + ' emptied because of --overwriteNanocorrect flag')
            else:
                try:
                    os.rmdir(nanocorrectDir) #remove it
                
                except OSError: #Cannot be removed, because it is not empty
                    raise RuntimeError('The nanocorrect directory is not empty! Is there a previous run present? --overwriteNanocorrect can be used to proceed', nanocorrectDir)
                #TODO: build a flag so this can be overwritten. Smart solution to store these files anyway
                logger.info(nanocorrectDir + ' existst but was empty, so proceed')
                os.makedirs(nanocorrectDir) #create it


    if not args.skipNanopolish: #Do a check for clean directory now instead of when demux is done
        ##Direcotry handeling
        nanopolishDir = getNanopolishDir(baseFileName)
        if not os.path.exists(nanopolishDir): #it it does not exist
            os.makedirs(nanopolishDir) #create it
            logger.info(nanopolishDir + ' did not exist, created it')
        else: #it exists
            logger.info(nanopolishDir + ' already exists')
            if args.overwriteNanopolish: #if we are sure to overwrite the existing nanopolish data?
                shutil.rmtree(nanopolishDir)
                os.makedirs(nanopolishDir) #create it
                logger.info(nanopolishDir + ' emptied because of --overwriteNanopolish flag')
            else:
                try:
                    os.rmdir(nanopolishDir) #remove it
                
                except OSError: #Cannot be removed, because it is not empty
                    raise RuntimeError('The nanopolish directory is not empty! Is there a previous run present? --overwriteNanopolish can be used to proceed', nanopolishDir)
                #TODO: build a flag so this can be overwritten. Smart solution to store these files anyway
                logger.info(nanopolishDir + ' existst but was empty, so proceed')
                os.makedirs(nanopolishDir) #create it
        if not os.path.exists(args.pathRawreads): #When we run nanopolish we need to have the raw reads defined
            raise IOError('the pathRawreads is not valid! This should point to your .fast5 files, see doc.')       
      
    if not args.skipCARD: #Do a check for clean directory now instead of when demux is done
        ##Direcotry handeling
        annotationDir = getAnnotationDir(baseFileName)
        if not os.path.exists(annotationDir): #it it does not exist
            os.makedirs(annotationDir) #create it
            logger.info(annotationDir + ' did not exist, created it')
        else: #it exists
            logger.info(annotationDir + ' already exists')
            if args.overwriteCARD: #if we are sure to overwrite the existing barcode output
                shutil.rmtree(annotationDir)
                os.makedirs(annotationDir) #create it
                logger.info(annotationDir + ' emptied because of --overwriteCARD flag')
            else:
                try:
                    os.rmdir(annotationDir) #remove it
                
                except OSError: #Cannot be removed, because it is not empty
                    raise RuntimeError('The CARD annotation directory is not empty! Is there a previous run present? --overwriteCARD can be used to proceed', annotationDir)
                #TODO: build a flag so this can be overwritten. Smart solution to store these files anyway
                logger.info(annotationDir + ' existst but was empty, so proceed')
                os.makedirs(annotationDir) #create it
      

    #####################################
    #Call all the relevant sub routines #
    #####################################

    
    if not args.skipDemux:
        deMux(baseFileName,args)
    else:
        logger.info('Skip the demux step because of --skipDemux')
    
    if not args.skipDemuxCollect: #a debug hook, skip demux but still use the .p file
        deMuxCollect(baseFileName,args)
    else:
        logger.info('Skip the demuxCollect step because of --skipDemuxCollect')
    
    if not args.skipNanocorrect:
        nanocorrect(baseFileName,args)
    else:
        logger.info('Skip the nanocorrect step because of --skipNanocorrect')

    if not args.skipNanopolish:
        nanopolish(baseFileName,args)
    else:
        logger.info('Skip the nanopolish step because of --skipNanopolish')
    

    if not args.skipCARD:
        annotateCARD(baseFileName,args)
    else:
        logger.info('Skip the CARD annotation step because of --skipCARD')
    
     
    
    logger.info('poreFUME done')
    logging.shutdown()

####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
# helper functions start. These are not directly called from main()
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################

def getBarcodeDir(baseFileName):
    return os.path.join('output','barcode',baseFileName)

def getNanocorrectDir(baseFileName):
    return os.path.join('output','nanocorrect',baseFileName)

def getNanocorrectDirABS(baseFileName):
    return os.path.join(os.getcwd(),'output','nanocorrect',baseFileName)

def getNanopolishDir(baseFileName):
    return os.path.join('output','nanopolish',baseFileName)

def getAnnotationDir(baseFileName):
    return os.path.join('output','annotation',baseFileName)


def spawnNanocorrect(thisInput):
    """
    The function spawns the nanocorrect pipeline, first the make and then the nanocorrect.py
    It takes a range of barcodes as argument that will be procced serially. This function is designed to be called using multiprocess so it can run in parallel.
    The round argument will not only decide the name of the interFixInput/Output but also if nanocorrect.py or nanocorrectC0.py is called
    Parameters:
    dictonary with:
        ranger = list of integers that represent that barcode to work on
        round =  round of nanocorrect we can do two round. so either  1 or 2
    
    """

    
    #unpack argument list
    thisRound = int(thisInput['round'])
    thisRange = thisInput['ranger']
    baseFileName = thisInput['baseFileName']
    pathNanocorrect = thisInput['pathNanocorrect']
    
    logger.debug('Spawned a process with range ' + str(thisRange) + ' for nanocorrect round ' + str(thisRound))
    
    if thisRound==1:
        interFixInput = 'afterBC'
        interFixOutput = 'afterNC1'
    elif thisRound==2:
        interFixInput = 'afterNC1'
        interFixOutput = 'afterNC2'
    else:
        raise ValueError('Invalid round number selected, only 1 or 2')
        
    
    for thisBarcode in thisRange: #Go through each barcode of the range that was passed
        
        logger.debug ('Running one-by-one nanocorrect on barcode ' + str(thisBarcode) + ' in round' +str(thisRound))
       
       #Make sure the input file can be used for error correction
        
        currentInputFile = os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode), ".".join([str(thisBarcode),baseFileName,interFixInput,'fasta'])  )       #makes output/nanocorrect/samplename/23/23.samplename.a
       
        
        if not os.path.isfile(currentInputFile): #Check the input fasta, ie the .afterBC. or afterNC1. really exists
            logger.warning(currentInputFile +  'does not exist. If this is a .afterNC1. file it can be nanocorrect did not generate output in the first round. Continue with next barcode')
            continue
        else:
            if len(list(SeqIO.parse(open(currentInputFile),"fasta"))) < 1: #Check how many records the start file has
                logger.warning(currentInputFile + ' does exist but has no records! If this is a .afterNC1. file it can be nanocorrect did not generate output in the first round. Continue with next barcode')
                continue
        
        
        
        ### RUN MAKE, this will run daligner under the hood to make the alignments (this is fast!)
        process = Popen(['make'
        ,'-f'
        ,os.path.join(pathNanocorrect,'nanocorrect-overlap.make')
        ,'INPUT='  + ".".join([str(thisBarcode),baseFileName,interFixInput,'fasta'])
        ,'NAME=' + str(thisBarcode)
    ], stdout=PIPE, stderr=PIPE, cwd= os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode)) ) #Set Current Working Directory to the barcode folder, as nanocorrect needs its own folder. This also means that the INPUT argument is relative to cwd and does not need a path
        
        stdout, stderr = process.communicate()
        process.wait() #wait till finished
        logger.debug(str(stdout))
        
        if stderr:
           logger.error('the nanocorrect MAKE file gave the following error while processing barcode: ' + str(thisBarcode) + ' : '+ str(stderr))
         
          
        
        logger.debug('Done with nanocorrect make for barcode ' +str(thisBarcode) + ' start nanocorrectX.py' + ' in round' +str(thisRound))
        
        
        ### RUN NANOCORRECT.py #Next step in the nanocorrect pipeline is to run nanocorrect.py itself. This calls poa on each alignment, which is a slow process
        if thisRound == 1:
            nanocorrectFilename = 'nanocorrect.py' #Has the min_coverage as 3
        elif thisRound == 2:
            nanocorrectFilename = 'nanocorrectC0.py' #Has the min_coverage as 1
        
        
        currentOutputFile = os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode), ".".join([str(thisBarcode),baseFileName,interFixOutput,'fasta'])  )       #makes
        
        correctFastaHandle = open(currentOutputFile, "wb") #Used to store the output
        process = Popen(['python'
        ,os.path.join(pathNanocorrect,nanocorrectFilename) #Executeable of nanocorret.py/nanocorrectC0.py
        , str(thisBarcode)
        ,'all'
    
    ], stdout=correctFastaHandle, stderr=PIPE, cwd= os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode)) ) #Set Current Working Directory to the barcode folder, as nanocorrect needs its own folder. This also means that the INPUT argument is relative to cwd and does not need a path
        
        stdout, stderr = process.communicate()
        process.wait() #wait till finished
        correctFastaHandle.close()
        logger.debug(str(stdout))
        
        if stderr:
           logger.error('the nanocorrect python script gave the following error while parsing barcode: ' + str(thisBarcode) + '.  '+ str(stderr))
        
        #Clean up the nanocorrect files
   
        os.remove(os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode), str(thisBarcode) + '.las'))
        os.remove(os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode),  str(thisBarcode) + '.pp.fasta.fai'))
        
        os.remove(os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode),  str(thisBarcode) + '.db'))
        os.remove(os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode),  str(thisBarcode) + '.pp.fasta'))
        
        os.remove(os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode),  '.' + str(thisBarcode) + '.idx'))
        os.remove(os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode),  '.' + str(thisBarcode) + '.dust.data'))
        os.remove(os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode),  '.' + str(thisBarcode) + '.dust.anno'))
        os.remove(os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode),  '.' + str(thisBarcode) + '.bps'))
        
        os.remove(os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode), 'HPCcommands.txt'))
        
        
        
        
        logger.debug('Done with nanocorrect for barcode ' + str(thisBarcode) + ' in round' +str(thisRound))
      
     



   

def getJobrange(totalBarcode,cores):
    """
    Distribute the barcodes over cores, can be used as range(start[i],end[i])
    Returns a (start,end)
    Paramters:
    totalreads = list of the barcodes to distribute, does not need to be consecutive ie. [1,3,9,12,48]
    cores = number of cores to split over ie 2
    
    """



    steps = int(math.ceil(len(totalBarcode)/float(cores)))
    start =  range(0,len(totalBarcode),steps)
    end = [x+steps for x in start]
    end[-1] = len(totalBarcode)
    
    return start,end
    
    #from http://biopython.org/wiki/Split_large_file
    #to avoid loading everything in memory
def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.
    
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
    
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch
            

def blastDatabase(queryFile,dbFile,args):
    """
    Returns a dataframe of a BLAST operation of the queryFile on the dbFile
    
    Parameters:
    queryFile = path to file that will be used as query in the BLAST process
    dbFile = path to file that will be used as database in the BLAST process
    """
    #Create a blast DB http://www.ncbi.nlm.nih.gov/books/NBK279688/

    
    assert type(queryFile) is StringType, "queryFile is not a string: %r" % queryFile
    assert type(dbFile) is StringType, "databaseFile is not a string: %r" % dbFile
    
    #Test inputs
    if not os.path.isfile(queryFile):
        raise IOError('queryFile does not excist',queryFile)
    
    if not os.path.isfile(dbFile):
        raise IOError('queryFile does not excist',dbFile)

    
    #Create dbase
    logger.info("Start builing blast database")
    process = Popen(['makeblastdb', '-in', str(dbFile), '-dbtype', 'nucl'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    process.wait() #wait till finished
    logger.info(str(stdout))
    if stderr:
        logger.error(str(stderr))
    
    #Search with blast
    logger.info("Start BLASTing for subsample %s",queryFile)
    process = Popen(['blastn','-db',str(dbFile),'-query',str(queryFile),'-max_hsps', '1', '-max_target_seqs', '1000', '-num_threads', str(args.cores), '-outfmt' ,'10','-out','blastn.tmp.output'], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    process.wait() #Wait till finished
    
    logger.info(str(stdout))
    if stderr:
        logger.error(str(stderr))
    logger.info("Finished BLASTing for subsample %s",queryFile)
    
    if not os.path.isfile('blastn.tmp.output'):
        raise RuntimeError('BLAST did not produce an alignment, one can implement an exception here. But for now stop')
    
    #Read in blast results
    dfBlast = None
    dfBlast = pd.read_csv('blastn.tmp.output',names=['qseqid' ,'sseqid' ,'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend' ,'sstart' ,'send','evalue', 'bitscore']) #load demuxed and corrected Minion dbase blasted by Sanger
    #TODO: make blastn baseFileName dependant!

    
    os.remove(str(dbFile) +'.nhr') #cleanup subsample file blastDB
    os.remove(str(dbFile) +'.nin') #cleanup subsample file blastDB
    os.remove(str(dbFile) +'.nsq') #cleanup subsample file blastDB
    
    os.remove('blastn.tmp.output') #cleanup blast output
    
    return dfBlast



def calcGeneLength(row):
    """
    The CARD database contains the position of the subject gene in the header name, this is originally stored in this script in
    id2 but split into a subjectGeneSTart and subjectGeneEnd. Based on this we can calculate the original gene length
    
    """
    
    return abs(int(row['subjectGeneStart'])- int(row['subjectGeneEnd']))




def calcCoverage(row):
    """
    We can calculate the coverage of the alignment by dividing the length of the alignemnt over the subjectGeneLength
    There is a glitch in the length of the header and the real length of the DNA in the card database, so limit of at 100%
    """
    #return (int(row['length'])/ float(row['subjectGeneLength']))*100
    
    return (int(row['length'])/ float(row['subjectGeneLength']))*100 if (int(row['length'])/ float(row['subjectGeneLength']))*100 < 100 else 100


def calcSegments(thisDF):
    """
    Calculates the most relevant hit for each segment on the read
    
    Paramters:
    thisDF = dataframe with the BLAST result for an INDIVIDUAL query. So don't pass the full BLAST table, but only from 1 qseqid, ie. thisDF[thisDF.qseqid == thisSeqid], where thisSeqid is the current seqid of interest
    """
    #from PIL import Image, ImageDraw #For viz purposes
    
    #im = Image.new('RGBA', (6000, 1000), (0, 0, 0, 0))  #initialize a debug drawing screen
    #draw = ImageDraw.Draw(im) 
    thisDF.reset_index(drop=True,inplace=True)

    dy = 50 #off set for debug drawing
    highscore = [] #Store the coordinates of the visited positions in a [start,end] format
    highindex = [] #Store the index number (=blast hit) of each visited position
    #for row in thisDF.itertuples(): #go through all the rows of the dataframe, not the most effient way! When this scales up write a new implemetation based on sorting

   
    for index,qstart,qend in zip(thisDF.index.values, thisDF.qstart, thisDF.qend): #this speeds up 5x compared to using iterrows(). However it is still slow. #TODO: make this vectorized if possible
        #print row
        if len(highscore) == 0: #First hit
            highscore.append([qstart,qend]) #Add position to the position list
            highindex.append(index)  #also keep track of the index
            #draw.line((row['qstart'],dy,row['qend'],dy), fill=(255,255,255)) #draw bright
        else:
            doesOverlap = False
            for thisScore in highscore: #Go through all the set positions
                if calcOverlap(thisScore[0],thisScore[1],qstart,qend) > 0: #If there is an overlap 
                    #print 'Overlapping, skip'
                    doesOverlap = True
                    break #Exit loop, no need to continue
            if doesOverlap == False: #We found a new segment
                highscore.append([qstart,qend]) #Add position to list
                #draw.line((row['qstart'],dy,row['qend'],dy), fill=(255,255,255))  #Draw bright
                highindex.append(index)

            else: #Old segment
                pass
                #draw.line((row['qstart'],dy,row['qend'],dy), fill=(55,55,55))     #Draw darker for debugging

        #print row['qstart'],row['qend']

        dy = dy + 2
    #im.save('out.png',"PNG")
    #im.show() #show a debug figure
 
    return thisDF.iloc[highindex] #Return a dataframe with the most relevant hit on each segment

def calcOverlap(a,b,c,d):
    """
    Returns the overlap between two segments (ab) vs (cd)
    
    Parameters:
    a = start position of segment AB
    b = end position of segment AB
    c = start position of segment CD
    d = end position of segment CD
    
    Can be in any direction and any order
    """
    return min([max([a,b]), max([c,d])]) - max([min([c,d]), min([a,b])]) #calculate overlap between two segments



def cmdExists(cmd):
    return any(
        os.access(os.path.join(path, cmd), os.X_OK)
        for path in os.environ["PATH"].split(os.pathsep)
    )


####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
# core functions start. These are directly called from main()
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################

def deMux(baseFileName,args):
    ###########################
    ##### Barcode Demux module
    ###########################
    """
    STEP 1a
    This basically finds the barcodes in the edges of a read using the smith waterman aligner.
    Results are retreived in a dataframe of which the barcodes are evaluated to pass a thershold (args.barcodeThreshold)
    """
    readcounter = 0
    
    scoreList = [2.7, -4.5, -4.7, -1.6]
    scoreList = [args.match,args.mismatch,args.gapopen,args.gapextend]
    shortHand = '_'.join([str(mli) for mli in scoreList])
    
    logger.info( 'Start with scoreList: ' + shortHand)
    logger.info('Split jobs over %s args.cores', args.cores)



    
    dictBarcode = {}
    
    for record in SeqIO.parse(args.fileBarcodes, "fasta"):
       
        
        seq = str(record.seq)
        
        if args.PacBioLegacyBarcode:
            ## This will add a 0 between the first 10 barcodes if it is entered like F_1 instead of F_01
            if record.id[-2] =="_": #from 1 to 9
                record.id = record.id[0:2] + '0' + record.id[-1]
        
        record.id = record.id[-2:]
        
        if record.id in dictBarcode: #lookup barcode in list
             dictBarcode[record.id].append(str(record.seq))
        else:
            dictBarcode.setdefault(record.id, [])
            dictBarcode[record.id].append(str(record.seq))
    
        
    
    #PB barcode is 5 padding + 16 barcodesc
    
    logger.info("Using barcode file: %s with %s barcodes",args.fileBarcodes,len(dictBarcode))
    logger.info("Using nanopore file: %s ",args.fileONTreads)
    logger.info("Runname: %s",baseFileName)
    logger.info("argument list: %s", args)


    
    record_iter = SeqIO.parse(open(args.fileONTreads),"fasta") #One generator for the file it self
    record_iter2 = SeqIO.parse(open(args.fileONTreads),"fasta") #And one generator to find the length
    
    inputLength = len(list(record_iter2))
    
    if inputLength<args.cores*2: #We run into trouble if there are more cores than sequences, since we get emtpy queues intitially. Since this is not a production run scenario we just stop the script here
        raise ValueError('There are less sequence records in the input file than there are cores defined. I guess this is a test run? Increase the amount of sequence in the input file, or lower the number of --cores')
    logger.debug("Each thread will have %s records to process",  int(math.ceil(float(inputLength)/args.cores)))
    #Split fasta files into batches so they can be processed parallel
    for i, batch in enumerate(batch_iterator(record_iter, int(math.ceil(float(inputLength)/args.cores)))):
        filename = baseFileName+ '.fasta.'+ str(i) +'.tmp'
        try:
            os.remove(filename)
        except OSError:
            pass
        handle = open(filename, "w")
        count = SeqIO.write(batch, handle, "fasta")
        handle.close()
    
    
    barcodeList = [] #Store barcode objects
    jobs = [] #Store multiprocessing objects
    
    #Run the barcoding in parallel
    for i in range(args.cores):
        
        thisFile = baseFileName+ '.fasta.'+ str(i) +'.tmp'
        thisInstance = Barcode(baseFileName,i)
        
        p = multiprocessing.Process(target=thisInstance.splitBarcode, args=(thisFile,dictBarcode, scoreList,args.barcodeEdge))
        
        jobs.append(p)
        p.start()
        barcodeList.append(thisInstance)
    
    for thisJob in jobs: #Wait till all jobs are done
        thisJob.join()


    
    
    #Finally collect the results from the barcoding in one dataframe
    dfCollector = pd.DataFrame()
    for i in range(args.cores):
        
        thisFilename = baseFileName + 'dfCollector.p.' + str(i) + '.tmp'
        
        dfthisCollector = pickle.load( open( thisFilename, "rb" ) )
        dfCollector = pd.concat([dfCollector,dfthisCollector])
        
        #cleanup
        os.remove(baseFileName+ '.fasta.'+ str(i) +'.tmp') #clenup temporary fasta file
        os.remove(thisFilename) #clean up temporary pickle file from processes

    
    count = str(len(dfCollector[dfCollector.wasMatch == 1].index))

    pickle.dump( dfCollector, open(os.path.join('output', baseFileName + ".afterBC.p"), "wb" ) )
    logger.info('Done demuxing, amount of double matches:' + str(len(dfCollector[dfCollector.wasMatch == 1].index)))

    
    logger.info('Start writing demuxed files')

        

def deMuxCollect(baseFileName,args):
    """
    Step 1b. Collect the results from the deMux() step
    """
    #Make sure the barcode directory is avialable
    
    barcodeDir  = getBarcodeDir(baseFileName)
    
    if not os.path.exists(barcodeDir): #it it does not exist
        os.makedirs(barcodeDir) #create it
        logger.info(barcodeDir + ' did not exist, created it')
    else: #it exists
        logger.info(barcodeDir + ' already exists')
        if args.overwriteDemux: #if we are sure to overwrite the existing barcode output
            shutil.rmtree(barcodeDir)
            os.makedirs(barcodeDir) #create it
            logger.info(barcodeDir + ' emptied because of --overwriteDemux flag')
        else:
            try:
                os.rmdir(barcodeDir) #remove it
            
            except OSError: #Cannot be removed, because it is not empty
                raise RuntimeError('The barcode directory is not empty! Is there a previous run present? --overwriteDemux can be used to proceed', barcodeDir)
            #TODO: build a flag so this can be overwritten. Smart solution to store these files anyway
            logger.info(barcodeDir + ' existst but was empty, so proceed')
            os.makedirs(barcodeDir) #create it
    
   
    
    dfCollector = pickle.load( open(os.path.join('output', baseFileName + ".afterBC.p"), "rb" ) )
    logger.info("%s amount of reads will be collected" ,dfCollector[(dfCollector[['scoreF','scoreR']].sum(axis=1)>args.barcodeThreshold)][['barcode']].stack().value_counts().sort_index().sum())
    counter = 0
    for record in SeqIO.parse(args.fileONTreads, "fasta"): #Walk through all the fasta sequences, this way requires a lot of writing IO, can be optimized by going per barcode
       
        logger.debug('Parsing: %s', record.description)
        
        if float(dfCollector[dfCollector.seqID == record.description][['scoreF','scoreR']].sum(axis=1)) > args.barcodeThreshold:
            logger.debug('Pass args.barcodeThreshold of ' + str(args.barcodeThreshold))
            if str(dfCollector[dfCollector.seqID == record.description]['direction'].item()) == 't':
                logger.debug('Strand is in template')
                start = int(dfCollector[dfCollector.seqID == record.description]['pos_F_end'])
                end = args.barcodeEdge-int(dfCollector[dfCollector.seqID == record.description]['pos_R_begin'])
                logger.debug('Forward cutoff:' +str(start) + ', end cutoff: -' + str(end))
                record.seq = record.seq[start:-end]
            elif str(dfCollector[dfCollector.seqID == record.description]['direction'].item()) == 'c':
                logger.debug('Strand is in complement')
                #For fun, but should write this in the manual as well we make the strand reverse complement.
                record.seq = bio.Seq.reverse_complement(record.seq)
                start = int(dfCollector[dfCollector.seqID == record.description]['pos_F_end'])
                end = args.barcodeEdge-int(dfCollector[dfCollector.seqID == record.description]['pos_R_begin'])
                logger.debug('Forward cutoff:' +str(start) + ', end cutoff: -' + str(end))
                record.seq = record.seq[start:-end]
            else:
                print str(dfCollector[dfCollector.seqID == record.description].direction.item())
                logger.error('Direction should be t or c')
                raise ValueError('Direction should be t or c')
            
            #open up directory
            
            #TODO: need to empty barcode files first!
            handle = open(  os.path.join(getBarcodeDir(baseFileName),str(int(dfCollector[dfCollector.seqID == record.description].barcode))+  '.' + baseFileName +'.afterBC.fasta'),"a")
            record.id = 'BC_' + str(int(dfCollector[dfCollector.seqID == record.description].barcode)) + '_' + record.id
            SeqIO.write(record,handle,"fasta")
            handle.close()
        else:
            handle = open(  os.path.join(getBarcodeDir(baseFileName), 'unknown.afterBC.fasta'),"a")
            SeqIO.write(record,handle,"fasta")
            handle.close()
            logger.debug('Lower (%s)then treshold, store as undefined' , float(dfCollector[dfCollector.seqID == record.description][['scoreF','scoreR']].sum(axis=1)))
        counter = counter + 1
        if counter%1000==0:
            logger.debug(counter)
            #break
            

def nanocorrect(baseFileName,args):
    """
    Step 2. Run two times error correction on the demuxed data
    
    """
    logger.info('Start Nanocorrect')
    
    
    #Load in the demuxed files
    
    barcodeFiles = next(os.walk(getBarcodeDir(baseFileName)))[2] #Find demuxed files in barode directory TODO:errorhandeling
    
    demuxedBarcodes = []
    logger.info( 'Work on the following files: ' + " ".join(barcodeFiles))
    
    if 'unknown.afterBC.fasta' in barcodeFiles: #Unknown barcodes will not be parsed
        barcodeFiles.remove('unknown.afterBC.fasta')
    
    for thisFilename in barcodeFiles:  #Go through each filename
        demuxedBarcodes.append(int(thisFilename.split('.')[0])) #Save the '8' of '8.poreCamp.2D.min500.afterBC.fasta'
    
    
    demuxedBarcodes =  sorted(demuxedBarcodes) #For fun, sort the list
    
   
    
    
    ### This returns a tuple of start and end. These numbers refer to the index of the demuxedBarcodes list. ie [(0, 11), (11, 22), (22, 33), (33, 43)]
    if len(demuxedBarcodes) == 0:
        raise ValueError('There are no barcodes to distribute. Did the deMux step result in any hits?')
        
    jobRange = zip(*getJobrange(demuxedBarcodes,args.cores))

    
    ###PREPARE THE DIRECTORY FOR NANOCORRECT
    for thisPosition in jobRange: #loop through each tuple of the barcode list
        
        thisRange = demuxedBarcodes[thisPosition[0]:thisPosition[1]]  #make a slice of the barcode list
        
        logger.debug('Ranges:')
        logger.debug (thisRange )
        
        for thisBarcode in thisRange:
            #Nanocorrect runs in an individual directory so need to create this
            os.makedirs(os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode)))
            
            
            #populate the individual barcode directory with the post barcode file
            #This will copy a file from output/barcode/mysample/32.mysample.afterBC.fasta to output/nanocorrect/mysample/32/32.mysample.afterBC.fasta
            copyfile(
            os.path.join(getBarcodeDir(baseFileName),".".join([str(thisBarcode),baseFileName,'afterBC.fasta'])) ,
            os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,'afterBC.fasta']) ))
           
    
    jobs = [] #Store multiprocessing jobs
    
    #### LOOP TO SPAWN IN PARALLEL NANOCORRECT round 1
    for thisPosition in jobRange: #loop through each tuple of the barcode
        
        thisRange = demuxedBarcodes[thisPosition[0]:thisPosition[1]]  #make a slice of the barcode list
        argList = {'round':1,
                    'ranger':thisRange,
                    'baseFileName':baseFileName,
                    'pathNanocorrect':args.pathNanocorrect}
        
        p = multiprocessing.Process(target=spawnNanocorrect, args=(argList,),) #denote list with ,!
        
        jobs.append(p)
        p.start()

    
    for thisJob in jobs: #Wait till all jobs are done
        thisJob.join()
    
    logger.info('Done with nanocorrect round 1! for all the barcodes ')
    
    jobs = []
    #### LOOP TO SPAWN IN PARALLEL NANOCORRECT round 2
    for thisPosition in jobRange: #loop through each tuple of the barcode
 
        
        thisRange = demuxedBarcodes[thisPosition[0]:thisPosition[1]]  #make a slice of the barcode list
        #Pack argument list
        argList = {'round':2,
                    'ranger':thisRange,
                    'baseFileName':baseFileName,
                     'pathNanocorrect':args.pathNanocorrect}
        
        p = multiprocessing.Process(target=spawnNanocorrect, args=(argList,),) #denote list with ,!
        
        jobs.append(p)
        p.start()

    
    for thisJob in jobs: #Wait till all jobs are done
        thisJob.join()
        
    
    logger.info('Done with nanocorrect round 2! for all the barcodes')
    
    
    
    ###Collect reads
    for thisDataset in ['afterBC','afterNC1','afterNC2']:
        
        storeRecords = []
        for thisBarcode in  demuxedBarcodes:
            try:
                handleInput = open(os.path.join(getNanocorrectDir(baseFileName),str(thisBarcode),'.'.join([str(thisBarcode), baseFileName, str(thisDataset), 'fasta'])), "r") #opens up output/nanocorrect/mysample/43/43.mysample.afterNC1.fasta
                        
            
            except IOError as errorCode:
                if errorCode.errno ==errno.ENOENT: #Source file doesnt exist, can happen if nanocorrect didnt yield any alignments
                    logger.warning('While collecting the results the followimg error occured:'+ str(errorCode))
                    continue
                else:
                    raise errorCode
            
            #append the barcode information in the readname
            for record in SeqIO.parse(handleInput, "fasta"):
                record.id = 'BC_' + str(thisBarcode) + '_' + record.description
                record.description = ''
                storeRecords.append(record)
            
            handleInput.close()
        
        handleOutput = open(os.path.join('output', baseFileName + '.'+thisDataset + '.fasta'), "w") #write to output/mysample.afterBC.fasta
        SeqIO.write(storeRecords, handleOutput, "fasta")
        handleOutput.close()
                
    
    logger.info('Done with nanocorrect module')
    
def nanopolish(baseFileName,args):
    """
    Step 3. Use Nanopolish to 'polish' the nanocorrect'ed reads.
    """

    #   Find barcode folders in nanocorrect
    #   popoulate nanopolish with BC and NC2
    #   make soft symlink to reads
    #   run nano polish 
    # collect reads?

    bwaDir = '/home/ubuntu/bwa/'
    nanopolishDir = '/home/ubuntu/nanopolish/'


    logger.info('Start Nanopolish')
    barcodes = next(os.walk(getNanocorrectDir(baseFileName)))[1] #Find the barcode directories (is index=1, for files, index=2) in the nanocorrect results

    
    logger.info( 'Work on the following files: ' + " ".join(barcodes))
    
    nanoBarcodes =  sorted([ int(x) for x in barcodes ]) #Convert strings to int and sort
    
   
    if len(nanoBarcodes) == 0:
        raise ValueError('There are no barcodes to distribute. Did the nanocorrect step result in any hits?')

    for thisBarcode in nanoBarcodes:
        logger.info("Start nanopolish loop on %s",thisBarcode)
        os.makedirs(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode))) #create a working directory

        tempFile2 = open("nanopolish.debug.2.txt","w")
        tempFile2.close()
        tempFile2 = open("nanopolish.debug.2.txt","a")
 

        tempFile = open("nanopolish.debug.1.txt","w")
        tempFile.close()
        tempFile = open("nanopolish.debug.1.txt","a")


        #copy in the models files, needed for R9, not R7.

        copyfile(nanopolishDir + 'etc/r9-models/nanopolish_models.fofn' ,
                    os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'nanopolish_models.fofn') )
        copyfile(nanopolishDir + 'etc/r9-models/template_median68pA.5mers.model' ,
                    os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'template_median68pA.5mers.model') )


        copyfile(nanopolishDir +'etc/r9-models/template_median68pA.model' ,
                    os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'template_median68pA.model') )


        #Nanopolish needs the raw event data. In order to get this the user supplies a pathRawreads to where the .fast5 files are located (currently only 1 raw read folder is supported). Next the first record is opened to find what paths is used to point the fast5 to by poretools. 
        for thisRead in SeqIO.parse(os.path.join(getNanocorrectDirABS(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,'afterBC.fasta'])),"fasta"): #Open the barcode split file, this should contain the info to the original read name
            fast5StorePlace =  os.path.split(thisRead.description.split(' ')[-1].rstrip())[0] #This will first split the header of the fasta record into chunks by ' '. The last chunk (-1) contains the path to the .fast5 file. To be sure remove the line end (rstrip). Then parse the filename and only take the 'head' and not the tail of the whole path. 
            break #we only need the first record, we assume the rest is the same
        
        #currenty a nested symlink is not supported
        if os.path.dirname(fast5StorePlace) == "":
            pass
        else:
            raise NotImplementedError('Currently poreFUME only supports 1-level deep directory pointers to the raw fast5 directory. ie. you ran poretools with poretools fasta --type 2D long/path/to/your/fast5/files/* instead run with poretools fasta --type 2D long/*')

        #instead of copying over all the raw read files (can be several hundereds of GB's) we make a symlink
        try:
            os.symlink(args.pathRawreads, os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),fast5StorePlace))
            
            #os.symlink('/home/ubuntu/extData/porecamp/pass/NB06/',os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'NB06'))
        except OSError as errorCode:
            if errorCode.errno==errno.EEXIST: #Symlink destination already exists
                logger.warning('Destination of symlink to raw read files already exist')
                os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),fast5StorePlace))
                #os.symlink('/home/ubuntu/extData/porecamp/pass/NB06/',os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'NB06'))
		os.symlink(args.pathRawreads, os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),fast5StorePlace))
            else:
                raise errorCode
              
                   
        #Make individual reads
        for thisRead in SeqIO.parse(os.path.join(getNanocorrectDirABS(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,'afterNC2.fasta'])), "fasta"):                    
            logger.info("Start extracting work of barcode: %s and read: %s\n",thisBarcode,thisRead.id)            
            readFile = open( os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta'])),"w")
            readFile.write(">" + str(thisRead.id) + '\n' )
            readFile.write( str(thisRead.seq) + '\n')
            readFile.close()

        # ../bwa/bwa index 37.poreCamp.2D.min500.OK.afterNC2.fasta
            runCmd = [os.path.join(bwaDir,'bwa'),
                                    'index',
                                    ".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta'])]
            p0 = Popen(runCmd,
                                    stdout=PIPE,stderr=tempFile,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))

            logger.info(" ".join(runCmd))
            print p0.communicate()
            p0.wait()

                   

            # ../bwa/bwa mem -x ont2d -t 4 37.poreCamp.2D.min500.OK.afterNC2.fasta 37.poreCamp.2D.min500.OK.afterBC.fasta | samtools view -Sb - | samtools sort -f - reads.sorted.bam
            runCmd = [os.path.join(bwaDir,'bwa'),
                                    'mem',
                                    '-x',
                                    'ont2d',
					'-t',
                                   str(args.cores) ,
                                    ".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta']),
                                    os.path.join(getNanocorrectDirABS(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,'afterBC.fasta']))]
            p1 = Popen(runCmd,
                                    stdout=PIPE,stderr=tempFile,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))
            logger.info(" ".join(runCmd))
            minAlign = int(len(thisRead.seq)*float(0.9))
            logger.info("Setting the minimum alignment to %s, of the total length %s",minAlign,len(thisRead.seq))
            runCmd = ['python', os.path.join(os.getcwd(),'piper.py'),
                                    '--minAlignment',
                                    str(minAlign)]
            p15 = Popen(runCmd,
                                    stdin=p1.stdout,stdout=PIPE,stderr=tempFile,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))
            logger.info(" ".join(runCmd))

            runCmd = ['samtools', 
                                    'view',
                                     '-Sb', 
                                     '-']
            p2 = Popen(runCmd, 
                                    stdin=p15.stdout, stdout=PIPE,stderr=tempFile,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))
            logger.info(" ".join(runCmd))


            runCmd = ['samtools', 
                                    'sort',
                                     
                                     '-o', 
                                     'reads.sorted.bam']
            p3 = Popen(runCmd,
                                    stdin=p2.stdout, stdout=PIPE,stderr=tempFile,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))
            logger.info(" ".join(runCmd))
            p3.wait()

            #samtools index reads.sorted.bam
            runCmd = ['samtools',
                                    'index',
                                    'reads.sorted.bam']
            p0 = Popen(runCmd,
                                    stdout=PIPE,stderr=tempFile,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))
            logger.info(" ".join(runCmd))
            print p0.communicate()
            p0.wait()

            proc = Popen(["samtools depth reads.sorted.bam | awk '{sum+=$3} END { print sum/NR}'"], stdout=PIPE, shell=True,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode))) ####!!!!! Running in shell !!!!!
            logger.info(" ".join(runCmd))
            try:
	        thisCoverage =  float(proc.communicate()[0])
            except:
                logger.error("Coverage depth of %s was not calculated succesfully, no reads mapped?. Continue with next read!",thisRead.id)
                thisCoverage = 0

 
            proc.wait()
            if thisCoverage<args.minCoverage:
                logger.warning("Coverage depth of %s is %s and thus lower than %s. Continue with next read",thisRead.id,thisCoverage,args.minCoverage)
                os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'reads.sorted.bam'))
                os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'reads.sorted.bam.bai'))
                os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta'])))
                os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta.bwt'])))
                os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta.pac'])))
                os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta.ann'])))
                os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta.amb'])))
                os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta.sa'])))
                continue
            else:
                logger.info("Coverage depth of %s is %s and thus higher than %s. Continue with this read",thisRead.id,thisCoverage,args.minCoverage)





            logger.info("Start nanopolish event align")
            #./nanopolish eventalign -t 4 --sam -r 37.poreCamp.2D.min500.OK.afterBC.fasta -b reads.sorted.bam -g 37.poreCamp.2D.min500.OK.afterNC2.fasta --models nanopolish_models.fofn | samtools view -Sb - | samtools sort -f - reads.eventalign.sorted.bam
            runCmd = [os.path.join(nanopolishDir,'nanopolish'),
                                    'eventalign' ,
                                    '-t',
                                    '1' ,
                                    '--sam' ,
                                    '-r',
                                    os.path.join(getNanocorrectDirABS(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,'afterBC.fasta'])),
                                    '-b' ,
                                    'reads.sorted.bam', #we can drop the abolulate path, as we use cwd=
                                    '-g' ,
                                    ".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta']),
                                    '--models',
                                    'nanopolish_models.fofn', #we can drop the abolulate path, as we use cwd=
                                    '-vvvv']
            pPolish = Popen(runCmd,
                                    stdout=PIPE,stderr=tempFile,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))


            logger.info(" ".join(runCmd))
            runCmd = ['samtools', 
                        'view',
                         '-Sb', 
                         '-']
            p2 = Popen(runCmd, 
                                    stdin=pPolish.stdout, stdout=PIPE,stderr=tempFile,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))

            logger.info(" ".join(runCmd))

            runCmd = ['samtools', 
                                    'sort',
                              
                                     '-o', 
                                     'reads.eventalign.sorted.bam']
            p3 = Popen(runCmd,
                                    stdin=p2.stdout, stdout=PIPE, stderr=tempFile, cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))
            logger.info(" ".join(runCmd))
            print p3.communicate()
            p3.wait()

            #samtools index reads.eventalign.sorted.bam
            runCmd = ['samtools',
                    'index',
                    'reads.eventalign.sorted.bam']
            p0 = Popen(runCmd,

                                    stdout=PIPE,stderr=tempFile,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))
            logger.info(" ".join(runCmd))
            print p0.communicate()
            p0.wait()

            logger.info("Done with nanopolish: event align\n")
            logger.info("Start with nanopolish parallelzation of barcode: %s and read: %s\n",thisBarcode,thisRead.id)            
            #python scripts/nanopolish_makerange.py 37.poreCamp.2D.min500.OK.afterNC2.fasta | parallel --results nanopolish.results -P 4 ./nanopolish variants --consensus polished.{1}.fa -w {1} -r 37.poreCamp.2D.min500.OK.afterBC.fasta -b reads.sorted.bam -g 37.poreCamp.2D.min500.OK.afterNC2.fasta -e reads.eventalign.sorted.bam -t 4 --min-candidate-frequency 0.1 --models nanopolish_models.fofn -vvvv
            runCmd = ['python' ,
                                   os.path.join(nanopolishDir,'scripts','nanopolish_makerange.py'),
                                   #'/Users/evand/Downloads/nanopolish/nanopolish/scripts/nanopolish_makerange.py',
                                    #os.path.join(getNanocorrectDirABS(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,'afterNC2.fasta']))
                                    ".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta'])
                                   ]
            p1 = Popen(runCmd,
                                    stdout=PIPE,stderr=tempFile,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))
            logger.info(" ".join(runCmd))
            runCmd = ['parallel',
                    '--gnu',
                                   '--results' ,
                                   'nanopolish.results',
                                   '-P',
                                   str(args.cores),
                                   os.path.join(nanopolishDir,'nanopolish'),
                                   'variants',
                                   '--consensus',
                                   'polished.{1}.' + thisRead.id +'.fa',
                                   '-w',
                                   '{1}',
                                   '-r',
                                   os.path.join(getNanocorrectDirABS(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,'afterBC.fasta'])),
                                   '-b',
                                   'reads.sorted.bam',
                                   '-g',
                                   ".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta']),
                                   '-e',
                                   'reads.eventalign.sorted.bam',
                                   '-t',
                                   '32',
                                   '--min-candidate-frequency',
                                   '0.1',
                                   '--models',
                                   'nanopolish_models.fofn',
                                   '']
            logger.info(" ".join(runCmd))
            try:
                p2 = Popen(runCmd, 
                                    stdin=p1.stdout, stdout=tempFile2,stderr=tempFile,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))
            except:
                logger.error("Somehow the Popen of nanopolish failed, this requires attention!")
                continue

            try:
                print p2.communicate()
            except:
                logger.error("Somehow the nanopolish command failed, this requires attention!")
                continue

            logger.info("Done with nanopolish parallelzation of barcode: %s and read: %s\n",thisBarcode,thisRead.id)            
            #if int(thisRead.id) == 1:
            #    logger.error("Reached programmed  end")
            #    break

            #Clean up all the intermidiate files
            os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'reads.eventalign.sorted.bam'))
            os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'reads.sorted.bam'))
            os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'reads.eventalign.sorted.bam.bai'))
            os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'reads.sorted.bam.bai'))
            os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta'])))
            os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta.bwt'])))
            os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta.pac'])))
            os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta.ann'])))
            os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta.amb'])))
            os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta.sa'])))
            os.remove(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),".".join([str(thisBarcode),baseFileName,thisRead.id,'fasta.fai'])))



        logger.info("Done with nanopolish parallelzation\n")
        logger.info("Printing logfile of the whole run now:\n")
        with open('tmpfile.txt') as f:
            for line in f:
                #print line
                logger.info(line)
        logger.info("Start collecting the polished data from all the barcodes\n")
        polishedFiles = []
        polishedOutput = open(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'allPolished.fasta'),"w")
        #print polishedFiles
        for thisFile in os.listdir(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode))):
    
            if thisFile.endswith(".fa") and thisFile.startswith("polished"):   
             polishedFiles.append(thisFile) #something is wrong with this intend! edit error?
          
           
        cmd = ['python' ,
                #'/Users/evand/Downloads/nanopolish/nanopolish/scripts/nanopolish_merge.py'
               os.path.join(nanopolishDir,'scripts','nanopolish_merge.py')
                               ]
        cmd.extend(polishedFiles)
        p1 = Popen(cmd,
                                stdout=polishedOutput,stderr=tempFile,cwd= os.path.join(getNanopolishDir(baseFileName),str(thisBarcode)))
        print p1.communicate()
            #python nanopolish_merge.py polished.*.fa > polished_genome.fa
            

    ###Collect polished data
    
    logger.info("Export nanopolish files from the allPolished to the .afterNP.")
    storeRecords = []
    for thisBarcode in  nanoBarcodes:
        try:
            handleInput = open(os.path.join(getNanopolishDir(baseFileName),str(thisBarcode),'allPolished.fasta'), "r") #opens up output/nanopolish/mysample/43/allPolished.fasta
                    
        
        except IOError as errorCode:
            if errorCode.errno ==errno.ENOENT: #Source file doesnt exist, can happen if nanopolish didnt yield any files
                logger.warning('While collecting the nanopolish results the followimg error occured:'+ str(errorCode))
                continue
            else:
                raise errorCode
        
        #append the barcode information in the readname
        for record in SeqIO.parse(handleInput, "fasta"):
            record.id = 'BC_' + str(thisBarcode) + '_' + record.description
            record.description = ''
            storeRecords.append(record)
        
        handleInput.close()

    handleOutput = open(os.path.join('output', baseFileName + '.afterNP.fasta'), "w") #write to output/mysample.afterBC.fasta
    SeqIO.write(storeRecords, handleOutput, "fasta")
    handleOutput.close()
    logger.info("nanopolish is done")

def annotateCARD(baseFileName,args):
    
    """
    Step 4. Use a BLAST against the CARD database to obtain an annotation
    
    """
    logger.info("Start CARD annotation")
    
    if not os.path.isfile(args.pathCARD):
        raise IOError('pathCARD does not exist! set with --pathCARD',args.pathCARD)

    assert type(baseFileName) is StringType, "baseFileName is not a string: %r" % baseFileName
    
    ###########
    ## Production loop
    ###########

    
    if args.annotateAll: #If al the datasets need to be annotated
        qpath = {
            'raw' :      os.path.join('inputData',  baseFileName + '.fasta'),
            'afterBC' :  os.path.join('output' , baseFileName + '.afterBC.fasta'),
            'afterNC1' : os.path.join('output', baseFileName + '.afterNC1.fasta'),
            'afterNC2' : os.path.join('output', baseFileName + '.afterNC2.fasta'),
  'afterNP' : os.path.join('output' ,baseFileName + '.afterNP.fasta')
        }
    else: #only annotate final set
        qpath = {
          'afterNP' : os.path.join('output' ,baseFileName + '.afterNP.fasta')
          ,'afterNC2' : os.path.join('output' ,baseFileName + '.afterNC2.fasta'),
    }
    
    #for thisQueryFile in qpath:
    for thisRunName, thisQueryFile in qpath.items():
       
        
        
        if not os.path.isfile(thisQueryFile):
            logger.warning('File to BLAST %s for CARD does not exist. Will not annotate and go to next.',thisQueryFile)
            continue
        
        recordIter = SeqIO.parse(open(thisQueryFile),"fasta") #And one generator to find the length
        
        inputLength = len(list(recordIter))
        
        if inputLength < 1:
            logger.warning('File to BLAST (%s) exists but contains no records! Will not annotate and go to next ',thisQueryFile)
            continue
        
        
        logger.info('Starting CARD analysis loop with %s stored at %s',thisRunName,thisQueryFile)
        dfBlast =  blastDatabase(thisQueryFile,args.pathCARD,args) #Run the blast command, also pass args so the function knows the cores
        
        #Parse the CARD database identifiers
        logger.info("Parse BLAST data")
        dfBlast[['gb','id1','id2','ARO','GeneName']]= dfBlast.sseqid.str.split('|',expand=True) #Split the CARD header out
        dfBlast[['subjectGeneStart','subjectGeneEnd']]= dfBlast.id2.str.split('-',expand=True) #Split the CARD gene length
        dfBlast = dfBlast.drop('id2',1) #remove id2 redundand
        dfBlast = dfBlast.drop('gb',1) #remove gb, does not contain information
        
        dfBlast['subjectGeneLength'] = dfBlast.apply(calcGeneLength, axis=1)
        
        dfBlast['coverageOfSubjectGene'] = dfBlast.apply(calcCoverage, axis=1)

        
        dfBlastSort = dfBlast.sort_values('bitscore',ascending=False) #sort by bitscore (in case it is not)
        
        dfBlastSort = dfBlastSort.reset_index(drop=True) #Reindex based on bitscore
        logger.info("Identified BLAST hits: " + str(dfBlastSort.shape[0]))
        
        logger.info("Start identifing relevant segments")

        
        #In order to calculate each the most relevant BLAST hits we go through each query sequence and calculate the relevant hits
        dfTotal = pd.DataFrame() #Initialze the final table
        
        for thisSeqid in dfBlastSort.qseqid.unique(): #go through each query sequence individually, here we could introduce multiprocessing
            
            dfTotal = dfTotal.append(calcSegments(dfBlastSort[dfBlastSort.qseqid == thisSeqid]), ignore_index=True) #Make a subset of the BLAST table containing only the current query seqience, next calculate the relevant blast hits and add them to the dfTotal table
        
        logger.info("Identified single hits: " + str(dfTotal.shape[0]))
        
        ###Save results in querySuffix + '.annotated.csv   
        dfTotal =  dfTotal.add_prefix('CARD:')
             
        dfTotal.to_csv(os.path.join(getAnnotationDir(baseFileName),os.path.splitext(os.path.basename(thisQueryFile))[0] + '.annotated.csv'))
        
        logger.info("Saved annotationed in: " + os.path.join(getAnnotationDir(baseFileName),os.path.splitext(os.path.basename(thisQueryFile))[0] + '.annotated.csv'))

        


if __name__ == "__main__":
    main()



