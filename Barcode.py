import logging
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
import Bio as bio

logger = logging.getLogger(__name__)    
 

class Barcode():
    
    thisunit = 'thisunit is at its initial value'
    def __init__(self,runname,thread=0):
        """
        Constructor
        Parameters:
        runname = name of the run, will be used to store tmp files
        thread = 0, will be used to store unique tmp files
        """
        logger.info('Barcode module initialized')
        self.runname = runname
        self.thread = thread


        
    def reverse(input):
        """
        Will give the reverse text of a string
        Parameters:
        input = input string
        """
        buildString = ''
        for x in input:
            buildString = x+buildString
        return buildString


 
    
    def splitBarcode(self,inputFasta,inputBarcode,scoreList,barcodeEdge):
        """
        Starts a loop through 1. the sequence file 2. the barcodes
        Parameters:
        inputFasta = path to fasta file containing nanopore sequenes
        inputBarcode = dict containing (asymetric) barcodes
        scoreList = list containing [match,mismatch,gapopen,gapextend] scores for alignment of barcode to sequence record
        barcodeEdge = distance from the read end to look for barcodes
        
        """
        logger.info('Start with %s in thread %s',inputFasta,self.thread)
        logger.info('Runname is %s',self.runname )
        listHighBarcodeScore = []
        listHighBarcodeID = []




        dfCollector = pd.DataFrame()

        for record in SeqIO.parse(inputFasta, "fasta"): #Walk through all the fasta sequences
        #for x in range(1,2):
            logger.debug('Parsing: %s', record.description)
            currentHigh = 0
            currentHighR = 0
           
            barcodeScore = {}
            barcodeScoreF = {}
            barcodeScoreR = {}


            df = pd.DataFrame()

            for thisBarcodeID,thisBarcodeSeq in inputBarcode.iteritems(): #Walk through the barcodes

                    alignmentsF = self.runPairwiseAlignment(str(record.seq)[0:barcodeEdge],thisBarcodeSeq[0],scoreList)
                    alignmentsR = self.runPairwiseAlignment(str(record.seq)[-barcodeEdge:],thisBarcodeSeq[1],scoreList)

                    alignmentsCF = self.runPairwiseAlignment(bio.Seq.reverse_complement(str(record.seq)[-barcodeEdge:]),thisBarcodeSeq[0],scoreList)
                    alignmentsCR = self.runPairwiseAlignment(bio.Seq.reverse_complement(str(record.seq)[0:barcodeEdge]),thisBarcodeSeq[1],scoreList)


        




                    #TODO: set barcode directly as index, prevents sizing up the dataframe every append
                    df = df.append({'barcode': thisBarcodeID, 
                                    'tF': alignmentsF[0][2], 
                                    'tR': alignmentsR[0][2], 
                                    'cF': alignmentsCF[0][2],
                                    'cR': alignmentsCR[0][2],
                                    'pos_tF_begin': alignmentsF[0][3],
                                    'pos_tF_end': alignmentsF[0][4],
                                    'pos_tR_begin': alignmentsR[0][3],
                                    'pos_tR_end': alignmentsR[0][4],
                                    'pos_cF_begin': alignmentsCF[0][3],
                                    'pos_cF_end': alignmentsCF[0][4],
                                    'pos_cR_begin': alignmentsCR[0][3],
                                    'pos_cR_end': alignmentsCR[0][4]


                                   }
                                   , ignore_index=True)

    

            df = df.set_index('barcode') #set the barcode as the index    
            #logger.info(df)
            hit = self.findHit(df)
            hit.update({'seqID': record.description})
            logger.debug(hit)

            dfCollector = dfCollector.append(hit,ignore_index=True)


            s = str(record.description)
        
        #print 'Saving pickle of thread %s' % thread
        dfCollector.to_pickle(self.runname + 'dfCollector.p.' + str(self.thread) + '.tmp')
        logger.info('Done with thread %s',self.thread)

        return True
        #return dfCollector

     
        
    def runPairwiseAlignment(self,inputSeq,barcodeSeq,scoring):
        """
        Does the actual calling of the pairwise2 algorithm
        """


     
        return pairwise2.align.localms(  str(inputSeq), str(barcodeSeq), *scoring)

      




    def getWinner(self,inputDf,column):
        """
        Get winner
        Returns list with
            barcode

        """
        hit = inputDf[[column]].sum(axis=1).sort_values(ascending=False).head(1).index[0]


      
        return [
            str(hit),
            float(inputDf[[column]].sum(axis=1).sort_values(ascending=False).head(1))
            ]



    def findHit(self,inputDf):
        """
        1. Checks if two barcodes are the top 1 score -> return this set
        2. If not, check which pair is highest scoring -> return this set

        Returns 
            barcode,
            direction, #complement or reverse
            scoreF, 
            scoreR,
            wasMatch, #case 1 (Above)
            rankF, #usefull in case 2 (in case one this is 1)
            rankR # "

        """




        if self.getWinner(inputDf,'tF')[0] == self.getWinner(inputDf,'tR')[0]:
            logger.debug("forward and reverse match in template strand")

            ##
            ## Check this piece, and then insert in complement
            ## 
            return {            'barcode': self.getWinner(inputDf,'tF')[0],
                                'direction': 't', 
                                'scoreF': self.getWinner(inputDf,'tF')[1], 
                                'scoreR': self.getWinner(inputDf,'tR')[1],
                                'wasMatch': True,
                                'rankF':1,
                                'rankR':1,
                                'pos_F_begin': int(inputDf.loc[self.getWinner(inputDf,'tF')[0]].pos_tF_begin),
                                'pos_F_end': int(inputDf.loc[self.getWinner(inputDf,'tF')[0]].pos_tF_end),
                                'pos_R_begin': int(inputDf.loc[self.getWinner(inputDf,'tR')[0]].pos_tR_begin),
                                'pos_R_end': int(inputDf.loc[self.getWinner(inputDf,'tR')[0]].pos_tR_end),

                    }


        if self.getWinner(inputDf,'cF')[0] == self.getWinner(inputDf,'cR')[0]:
            logger.debug("forward and reverse match in compelement strand")
            ##
            ## Check this piece
            ## 
            return {            'barcode': self.getWinner(inputDf,'cF')[0],
                                'direction': 'c', 
                                'scoreF': self.getWinner(inputDf,'cF')[1], 
                                'scoreR': self.getWinner(inputDf,'cR')[1],
                                'wasMatch': True,
                                'rankF':1,
                                'rankR':1,
                                'pos_F_begin': int(inputDf.loc[self.getWinner(inputDf,'cF')[0]].pos_cF_begin),
                                'pos_F_end': int(inputDf.loc[self.getWinner(inputDf,'cF')[0]].pos_cF_end),
                                'pos_R_begin': int(inputDf.loc[self.getWinner(inputDf,'cR')[0]].pos_cR_begin),
                                'pos_R_end': int(inputDf.loc[self.getWinner(inputDf,'cR')[0]].pos_cR_end)
                   }
        else:
            logger.debug("no double forward and reverse match found")





        topT = inputDf[['tF','tR']].sum(axis=1).sort_values(ascending=False).head(1)
        topC = inputDf[['cF','cR']].sum(axis=1).sort_values(ascending=False).head(1)

        if float(topT) >= float(topC): #Template score is higher than complement score
            logger.debug( 'Template strand has highest scoring duo')
            winner = inputDf.ix[inputDf[['tF','tR']].sum(axis=1).sort_values(ascending=False).index].head(1)
            #print df[['tF','tR']].ix[df[['tF','tR']].sum(axis=1).sort_values(ascending=False).index]
            rankWinner = inputDf.rank(axis=0,            #over the rows
                                 numeric_only=True, #only values, sanaty check should render this not needed
                                 ascending =False   #high to low
                                ).ix[winner.index.tolist()[0]]
            #print df.rank(axis=0,            #over the rows
            #                     numeric_only=True, #only values, sanaty check should render this not needed
            #                     ascending =False   #high to low
            #                    )
            #print 'Winner barcode:',
            #print winner.index.tolist()[0]
            #print 'Rank dataframe:'
            #print df.rank(axis=0,numeric_only=True,ascending=False )
            #print 'Rank of the winner:'
            #print rankWinner
            return {    'barcode': winner.index.tolist()[0],
                        'direction': 't', 
                        'scoreF': float(winner['tF']), 
                        'scoreR': float(winner['tR']),
                        'wasMatch': False,
                        'rankF':rankWinner['tF'], 
                        'rankR':rankWinner['tR'],
                        'pos_F_begin':int(winner.pos_tF_begin),
                        'pos_F_end': int(winner.pos_tF_end),
                        'pos_R_begin': int(winner.pos_tR_begin),
                        'pos_R_end': int(winner.pos_tR_end)
           }

        else:
            logger.debug( 'Complement strand has highest scoring duo')
            winner = inputDf.ix[inputDf[['cF','cR']].sum(axis=1).sort_values(ascending=False).index].head(1)
            rankWinner = inputDf.rank(axis=0,            #over the rows
                                 numeric_only=True, #only values, sanaty check should render this not needed
                                 ascending =False   #high to low
                                ).ix[winner.index.tolist()[0]]




            return {    'barcode': winner.index.tolist()[0],
                        'direction': 'c', 
                        'scoreF': float(winner['cF']), 
                        'scoreR': float(winner['cR']),
                        'wasMatch': False,
                        'rankF':rankWinner['cF'], #todo
                        'rankR':rankWinner['cR'],
                        'pos_F_begin': int(winner.pos_cF_begin),
                        'pos_F_end': int(winner.pos_cF_end),
                        'pos_R_begin': int(winner.pos_cR_begin),
                        'pos_R_end': int(winner.pos_cR_end)
           }


