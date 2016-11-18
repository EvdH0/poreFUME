# Installing poreFUME

The easiest way is to run:

1.  ```install.sh``` which takes care of [daligner](https://github.com/thegenemyers/DALIGNER), [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB), [POA](http://sourceforge.net/projects/poamsa/), [nanocorrect](https://github.com/jts/nanocorrect) and BLAST+
2.  register the path variables by running ```. env.sh``` in the current shell. You can export the path in your ```~/.profile```. 

You will probably need a few dependencies as listed below.

## Requirements

poreFUME requires Python 2.7 or newer. Furthermore the following packages can be installed accordingly
```
 apt-get install git
 apt-get install make
 
 apt-get install python-nose
 apt-get install python-pip
 apt-get install python-pip
 apt-get install python-dev

 apt-get install bioperl

 apt-get install build-essential
 apt-get install zlib1g-dev

 apt-get install libncurses5-dev #needed for samtools
  
 pip install biopython
 pip install numpy
 pip install pandas
 pip install cython
 pip install pysam
```

##Testing 
To test the working of poreFUME you can run ```nosetests -v``` which should output something like
```
ubuntu@testhost:~/poreFUME$ nosetests -v
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
