poreFUME
===============

Demultiplex, correct and annotate antibiotic resistance genes in nanopore data




# Usage

```
usage: poreFUME.py [-h] [--PacBioLegacyBarcode] [--verbose] [--overwriteDemux]
                   [--overwriteNanocorrect] [--overwriteCARD] [--skipDemux]
                   [--skipDemuxCollect] [--skipNanocorrect] [--skipCARD]
                   [--match [MATCH]] [--mismatch [MISMATCH]]
                   [--gapopen [GAPOPEN]] [--gapextend [GAPEXTEND]]
                   [--cores [CORES]] [--barcodeThreshold [BARCODETHRESHOLD]]
                   [--barcodeEdge [BARCODEEDGE]]
                   [--pathNanocorrect [PATHNANOCORRECT]]
                   [--pathCARD [PATHCARD]] [--annotateAll]
                   fileONTreads fileBarcodes

positional arguments:
  fileONTreads          path to FASTA where the (2D) nanopore reads are stored
  fileBarcodes          path to FASTA where the barcodes are stored, format
                        should be ie F_34 for forward and R_34 for reverse
                        barcode

optional arguments:
  -h, --help            show this help message and exit
  --PacBioLegacyBarcode
                        the pacbio_barcodes_paired.fasta file has first digist
                        as 4 instead of 04, turning this option on will fix
                        this
  --verbose             switch the logging from INFO to DEBUG
  --overwriteDemux      overwrite results in the output/barcode/runid
                        directory if they exist
  --overwriteNanocorrect
                        overwrite the results in the output/nanocorrect/runid
                        directory if the exist
  --overwriteCARD       overwrite the results in the output/annotation/runid
                        directory if the exist
  --skipDemux           Skip the barcode demux step and proceed with
                        nanocorrect, cannot be used with overwrite. Assumes
                        the output/barcode/ and output/ directory are
                        populated accordingly
  --skipDemuxCollect    will skip the demux it self and go to colletion based
                        on the pickle
  --skipNanocorrect     Skip the nanocorrect step.
  --skipCARD            Skip the CARD annotation
  --match [MATCH]       Score for match in alignment (default: 2.7)
  --mismatch [MISMATCH]
                        Score for mis-match in alignment (default: -4.5)
  --gapopen [GAPOPEN]   Score for gap-open in alignment (default: -4.7)
  --gapextend [GAPEXTEND]
                        Score for gap-extend in alignment (default: -1.6)
  --cores [CORES]       Amount of args.cores to use for multiprocessing
                        (default: 1)
  --barcodeThreshold [BARCODETHRESHOLD]
                        Minimum score for a barcode pair to pass (default: 58)
  --barcodeEdge [BARCODEEDGE]
                        Maximum amount of bp from the edge of a read to look
                        for a barcode. (default: 60)
  --pathNanocorrect [PATHNANOCORRECT]
                        Set the path to the nanocorrect files (default:
                        /Users/evand/Downloads/testnanocorrect/nanocorrect/)
  --pathCARD [PATHCARD]
                        Set the path to CARD fasta file (default:
                        inputData/n.fasta.protein.homolog.fasta)
  --annotateAll         By default only the final (demuxed and two times
                        corrected) dataset is annotated, however by turning on
                        this option all the files, raw, after demux, after 1st
                        round of correction, after 2nd round of correction are
                        annotated. This obviously takes longer.
```


### Example

```poreFUME.py inputData/2DnanoporeData.fasta inputData/barcodes.fasta --PacBioLegacyBarcode --barcodeThreshold 50 --annotateAll --verbose --cores 8```

## Output :

Folders ```output/barcode/mySample/```, ```output/nanocorrect/mySample/```, ```output/annotation/yourSample/``` ```and output/``` contain the output files. 

The pipeline consists of 3 steps:
1.  barcode demultiplexing
2.  error correction using nanocorrect
3.  annotation of the reads using the CARD datbase

Each step can be skipped using the relevant ```--skipXXX``` parameter. For example when only barcodes need to extracted ```--skipNanocorrect``` and ```--skipCARD``` can be used.

### Step 1. demultiplexing of barcodes

When ```--skipDemux``` is not set (_default_), ```output/mySample.afterBC.p``` will be created which contains a pickeled pandas dataframe with the barcode score for each read. Relevant parameters are ```--match```,```--mismatch```,```--gapopen``` and ```--gapextend``` which can be used to adjust the score function of the barcode alignment.  

When ```--skipDemuxCollect``` is not set (_default_): ```output/mySample.afterBC.fasta``` will be created based on the ```output/mySample.afterBC.p``` file. The reads in the ```output/mySample.afterBC.fasta``` are identified by the ```FASTA header >BC_{barcodeID}_{originalreadname}``` ie. ```>BC_39_nanporeEcoliGenomeReadHash-3a43-4j34...```. Furthermore ```output/barcode/mySample/``` will contain a FASTA file for each individual barcode, again with the barcode in the FASTA header. Reads on which barcodes could not be accurately  determined given a ```--barcodeThreshold``` are placed in unknown.fasta.   
```--barcodeEdge``` is used to determine how far barcodes are searched in the read, increasing this value linearly scales with the run time of the step.  
_Note: when poreFUME already find data in ```output/barcode/mySample/``` it will terminate, this can be overruled by passing the ```--overwriteDemux``` flag. This will remove all the existing data in ```output/barcode/mySample/```_.

### Step 2. error correction of the demultiplexed reads using nanocorrect
When ```--skipNanocorrect``` is not set (_default_) poreFUME will invoke [nanocorrect](https://github.com/jts/nanocorrect) to error correct the demultiplexed reads. Since nanocorrect needs its own directory to run in when ran in parallel, it will create ```/output/nanocorrect/mySample/{barcodeID}/``` directories. Inside [DALIGNER](https://github.com/thegenemyers/DALIGNER) and [poa](https://sourceforge.net/projects/poamsa/) will be run as called by nanocorrect. ```--pathNanocorrect``` can be set to point to the nanopore package.  
_Note: when poreFUME already find data in ```output/nanocorrect/mySample/``` it will terminate, this can be overruled by passing the ```--overwriteNanocorrect``` flag. This will remove all the existing data in ```output/nanocorrect/mySample/```_.

### Step 3. annotation using CARD
The final step is annotation of the data using the [CARD database](https://card.mcmaster.ca/) when ```--skipCARD``` is not set (_default_). This part will look for ```output/mySample.afterNC2.fasta``` and annotate the file against the CARD database.  ```--pathCARD``` is used to point to the nucleotide CARD database. When ```--annotateAll``` is set, also ```inputFiles/mySample.fasta``` , ```output/mySample.afterBC.fasta``` and ```output/mySample.afterNC1.fasta``` will be annotated. The output is a CSV file in ```output/annotation/mySample/mySample.AFTERNC2.annotation.csv``` containing the readname and the relevant CARD information.   
_Note: with ```--skipCARD``` the output in ```output/annotation/mySample/``` will be overwritten_.

### Parallelization
With the ```--cores``` flag the following processes can be parallelized: 
* Smith-Waterman algorithm to detect barcodes
* nanocorrect on multiple barcodes
* BLAST to annotate with the CARD database. 

## Testing
To test the working of poreFUME you can run ```nosetests -v``` which should output something like
```
ubuntu@ip-172-3:~/poreFUME$ nosetests -v
test this ... ok
testCARDavialable (test.TestCARD) ... ok
testInputavialable (test.TestCARD) ... ok
testSegments (test.TestCARD) ... ok
testBLASTDATABASE (test.TestDependencies) ... ok
testBLASTN (test.TestDependencies) ... ok
testDBdust (test.TestDependencies) ... ok
testDBsplit (test.TestDependencies) ... ok
testF2DB (test.TestDependencies) ... ok
testLAcat (test.TestDependencies) ... ok
testPOA (test.TestDependencies) ... ok
job ranger returns index of begin and end of job range. ... ok
testOverlap (test.TestFunctions) ... ok

----------------------------------------------------------------------
Ran 13 tests in 0.272s

OK
```

## Requirements
poreFUME makes use of the [CARD database](https://card.mcmaster.ca/). So when using please cite [McArthur et al. 2013. The Comprehensive Antibiotic Resistance Database. Antimicrobial Agents and Chemotherapy, 57, 3348-3357.](http://www.ncbi.nlm.nih.gov/pubmed/23650175). Furthermore [Nanocorrect](https://github.com/jts/nanocorrect) is used, which can be cited by [Loman NJ, Quick J, Simpson JT: A complete bacterial genome assembled de novo using only nanopore sequencing data. Nat Methods 2015, 12:733–735.](http://www.nature.com/nmeth/journal/v12/n8/abs/nmeth.3444.html)
