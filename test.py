
import poreFUME
import unittest
import os



from argparse import Namespace


import pandas as pd
from pandas.util.testing import assert_frame_equal   

class TestCARD(unittest.TestCase):

 
    def testCARDavialable(self):
        assert os.path.exists(os.path.join('test','data','n.fasta.protein.homolog.fasta'))
    
    def testInputavialable(self):
        assert os.path.exists(os.path.join('test','data','CblA_barcode48.fasta'))
           
    
    def testBlast(self):
        """
        test this
        """
        args = Namespace(cores=1)
        returnDF = poreFUME.blastDatabase(os.path.join('test','data','CblA_barcode48.fasta'),os.path.join('test','data','n.fasta.protein.homolog.fasta'),args)
        self.assertEqual(returnDF.shape[0],1)
        self.assertEqual(str(returnDF['sseqid'][0]),'gb|GQ343019|132-1023|ARO:3002999|CblA-1')
        
    
    def testSegments(self):
        inputDF = pd.read_pickle(os.path.join('test','data','consensus9.input.df.p'))
        outputDF = pd.read_pickle(os.path.join('test','data','consensus9.output.df.p'))
  
        assert_frame_equal(poreFUME.calcSegments(inputDF),outputDF) #Test dataframes are 

class TestFunctions(unittest.TestCase):

    def testJobrange(self):
        """
        job ranger returns index of begin and end of job range. 
        """
        self.assertEqual(poreFUME.getJobrange([01,11,88,99,1987,2011,3044,6789],3),([0, 3, 6], [3, 6, 8]))
    
    def testOverlap(self):
        self.assertEqual(poreFUME.calcOverlap(10,25,20,60),5)
        self.assertEqual(poreFUME.calcOverlap(10,25,30,60),-5)
        
        

class TestDependencies(unittest.TestCase):

    def testPOA(self): 
        self.assertTrue(poreFUME.cmdExists('poa'))
        
    def testBLASTN(self):
        self.assertTrue(poreFUME.cmdExists('blastn'))
    
    def testBLASTDATABASE(self):
        self.assertTrue(poreFUME.cmdExists('makeblastdb'))
    
        
    def testF2DB(self):
        self.assertTrue(poreFUME.cmdExists('fasta2DB'))    
    
    def testDBsplit(self):
        self.assertTrue(poreFUME.cmdExists('DBsplit'))  
        
    def testDBdust(self):
        self.assertTrue(poreFUME.cmdExists('DBdust'))
        
    def testLAcat(self):
        self.assertTrue(poreFUME.cmdExists('LAcat'))
if __name__ == '__main__':
    
    
 
    unittest.main()
