#! /usr/bin/env python
'''
cluster.py
Uses pamk from clump.py to determine clustering for mutation positions in a given gene
'''

from optparse import OptionParser
#from clump import pamk
import rpy2.robjects as robjects
from Bio.PDB.PDBParser import PDBParser
import math
import pandas as pd
import sys, os
import numpy as np
import rpy2.robjects.numpy2ri




#import random

r=robjects.r
r.library('fpc')

rpy2.robjects.numpy2ri.activate()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
class CmdOpts(object):
    usage="""%prog [options] -f -p


    """

    def __init__(self):
        parser = OptionParser(usage=CmdOpts.usage)
        parser.add_option("-f", "--fname", dest="filename",
                          help="""Input annotation file""")
        parser.add_option("-a","--aname",dest="affreq",
                          help="""Input allele frequency cutoff""",default=1)
        parser.add_option("-p","--pname",dest='proteinlengths',
                          help="""Input file of NP and protein lengths""")
        parser.add_option('-c','--cname',dest='controlperms',
                          help='''Input the controls that will be used for permutation testing''')
        parser.add_option('-z','--zname',dest='permutations',
                          help='''Input the number of permutations you want to perform on each gene''')

        parser.add_option('-s','--sname',dest='structurepath',
                          help='''Input the file directory of the pdb structures you want to use''')

        parser.add_option("-m","--mname",dest='nmutcutoff', default=5,
                          help="""Number of mutations needed in a gene to be considered""")
        parser.add_option('-n',action="store_true",dest='normalize',
                          help="""Normalize protein lengths""")
        parser.add_option('-t',action='store_true',dest='titles',help="""Output Column Titles""")
        
        (opts, args) = parser.parse_args()

        if not opts.filename and not opts.affreq and not opts.proteinlengths:
            parser.error("Incorrect input args")

        self.__dict__.update(opts.__dict__)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
class Parser():
    '''
    Parses the AD_final and AR_final files
    '''
    def FillLen(self,Proteinfile):
        lendb=dict()
        for line in open(Proteinfile):
            lineinfo=line.strip().split('\t')
            lendb[lineinfo[0]]=int(lineinfo[1])
        return lendb
    
    def get_coord(self, pos, protein_id):
        '''
        Returns the 3d coordinates for a mutation in a given protein
        '''

        return self.db[protein_id]

    def fill(self, testfile,afcutoff):
         # get pdb files for each protein
        pdbparser = PDBParser()

        db=dict()
        domain=dict()
        for line in open(testfile):
            
            lineinfo = line.strip().split('\t')
            key=(lineinfo[0],lineinfo[1],lineinfo[2])
            
            if len(lineinfo)==8:
                allelefreq=0
            elif lineinfo[8]=="NA":
                allelefreq=0
            else:
                allelefreq=float(lineinfo[8])
            if allelefreq <= afcutoff:
                pos = int(lineinfo[3])
                protID = lineinfo[1]

                if protID not in self.pbds:
                    try:
                        self.pbds[protID] = pdbparser.get_structure(protID, self.structdir + protID + '.pdb')[0]['A']
                    except FileNotFoundError:
                        continue
            

                ########
                #structure = pdbparser.get_structure(protID, "/Users/elaine/Library/CloudStorage/OneDrive-JohnsHopkins/CLUMP-master/structures/" + protID + ".pdb")
                structure = self.pbds[protID]
                # get the coordinates for the mutation
                #print(protID, pos)
               	coord = np.array(structure[pos]['CA'].get_coord())
           

                ########
                

                if self.normalize:
                    plen = self.Proteinlen[protID]
                    pos=float(pos)/float(plen)
                

                if key in db.keys():
                    db[key] = np.vstack((db[key], [coord]))
                    
                    if len(lineinfo)==9:
                        if lineinfo[9]!="NA" and lineinfo[9]!="NONE":
                            domain[protID]=domain[protID]+1
                else:
                    db[key]=[coord]
                    if len(lineinfo)==9:
                        if lineinfo[9]!="NA" and lineinfo[9]!="NONE":
                            domain[protID]=1
                    else:
                        domain[protID]=0

        return db,domain

    def __init__(self, inputfile,lengthfile,af,ncutoff,normalize,structdir):
        

        self.Proteinlen = self.FillLen(lengthfile)
        self.normalize=normalize
        self.NMUTCUTOFF=ncutoff
        self.structdir=structdir
        #Initializing dictionary of key : (pos, coord)
        self.db = dict()
        #domain db
        self.domain=dict()
        #Initializing medoids
        self.medoids = dict()
        #Initializing medoids dict for easy calling (keyed by NPid)
        self.NP2medoids = dict()
        #Initializing clump scores
        self.clumps = dict()

        #Initializing pdb structures
        self.pbds = dict()

        #Now reading in file:
        self.db,self.domain = self.fill(inputfile,af)
        #Now that the dictionary is filled, applying clustering to the mutation positions
        self.medoids = dict()
        #Grabbing the medoids for this particular key
        for key in self.db:
            sys.stderr.write('Now clustering: ' + str(key[0]) + '\n')
            muts = self.db[key]
            
            if len(muts) < int(self.NMUTCUTOFF): continue
            medoids = self.pamk(self.db[key])
            self.medoids[key] = medoids
            protID = key[1]
            self.NP2medoids[protID] = medoids

        #Applying clump
        self.clump()


    def __call__(self):
        '''
        Prints out the clump score for each record
        '''
        for record in self.clumps:
            print('\t'.join(record) + '\t' + str(self.clumps[record]))

    def getmedoids(self, protID):
        try:
            return self.NP2medoids[protID]
        except KeyError:
            return None
        

    def pamk(self,data):
        '''
        Wrapper around the fpc library's pamk function '''

        # now data is a 3*#mutations list 
        #data = [float(x) for x in data]

        '''
        if len(data) <= 2:
            mean = sum(data)/len(data)
            return [mean]
        '''
        #print(data)
        r = robjects.r
        nr,nc = data.shape
        m = r.matrix(data, nrow=nr, ncol=nc)
        r.assign("data", m)
        #m = r.matrix(data, nrow = len(data), ncol = 3)
        #m = robjects.r.matrix(robjects.FloatVector(data), byrow = True, ncol=3)
        if len(data) < 4:
            krange = [2]
        else:
            kmax = int(round(len(data)/2.0))
            krange = range(2,min(kmax+1, 15))
        try:
            if len(data) > 1000:
                reply = r.pamk(m, robjects.IntVector(krange), usepam = False)
                pamobject = reply[0]
                medoids = list(pamobject[1])
            else:
                reply = r.pamk(m, robjects.IntVector(krange))
                pamobject = reply
                pamobject = pamobject[0]  # mediod as r matrix
                medoids = np.array(pamobject[0])
            return medoids
        except:
            raise
            return None


    def clump(self):
        '''
        Applies the CLUMP score for all records
        This CLUMP score is the average of the log distance of a mutation to its nearest medoid
        '''
        sys.stderr.write('Now clumping\n')
        
        for record in self.medoids:
            
            distances = []
            medoids = self.medoids[record]
            muts = self.db[record]
            #muts = [coords[i:i+3] for i in np.arange(0, len(coords), 3)]
            for mutpos in muts:
            
                dist = math.log(min(np.sqrt(np.sum((medoids - mutpos)**2, axis = 1))) + 1)
                #dist =  math.log(min([np.sqrt(np.sum((mutpos-x)**2)) for x in medoids]) + 1)
                #dist = math.log( min( [abs(mutpos - x) for x in medoids] ) + 1)
                distances.append(dist)
            #Now taking the average
            clump = float(sum(distances)) / float(len(muts))
            self.clumps[record] = clump


    def permparse(self,filehandle,afcutoff):
        '''parse file to get db for permutations'''
        '''parse file to get db for permutations'''
        db=dict()
        for line in open(filehandle):
            lineinfo = line.strip().split('\t')
            key=(lineinfo[0],lineinfo[1],lineinfo[2])
            allelefreq=float(lineinfo[8])
            if allelefreq<=afcutoff:
                protID = lineinfo[1]
                if protID in db.keys():
                    db[protID]=db[protID]+1
                else:
                    db[protID]=1
        return db



    def controlclump(self,db):
        '''perform clump on the control population'''

        clumps=dict()
        for key in db:
            medoids=self.pamk(db[key])
            distances = []
            muts = db[key]
            #muts = [coords[i:i+3] for i in np.arange(0, len(coords), 3)]
            for mutpos in muts:
                dist = math.log(min(np.sqrt(np.sum((medoids - mutpos)**2, axis = 1))) + 1)
                #dist = math.log( min( [abs(mutpos - x) for x in medoids] ) + 1)
                distances.append(dist)
            #Now taking the average                                                                                             
            clump = float(sum(distances)) / float(len(muts))
            clumps[key] = clump
        return clumps


    def permutedb(self,db,numberofperm):
        '''do permutations'''
        pdbparser = PDBParser()
        clumppermutations=dict()
        for key in db:
            #print(key)
            protID = key
            #structure = pdbparser.get_structure(protID, "/Users/elaine/Library/CloudStorage/OneDrive-JohnsHopkins/CLUMP-master/structures/" + protID + ".pdb")
            structure = self.pbds[protID]

            if db[key]>=int(self.NMUTCUTOFF):##this was added to correct when controls have less than ncutoffs
                clumppermutations[key]=[]
            
                for ii in range(0,numberofperm):
                    #print(self.Proteinlen, db)
                    try: 
                        # get list of mutated positions                                      
                        permsample=np.random.choice(range(1, self.Proteinlen[key] + 1),size=db[key],replace=True)
                        # get the coordinates for the mutation
                        perm_coords = np.empty((0, 3))
                        for pos in permsample:
                            #print(pos)
                            perm_coords = np.vstack((perm_coords, [list(structure[int(pos)]['CA'].get_coord())]))
           
                    
                        #permsample = perm_coords 
                        permmedoids=self.pamk(perm_coords)

                        #permsample = [coords[i:i+3] for i in np.arange(0, len(coords), 3)]
                        distances=[]
                        for perm in perm_coords:
                            dist = math.log(min(np.sqrt(np.sum((permmedoids - perm)**2, axis = 1))) + 1)
                            #dist =  math.log(min([np.sqrt(np.sum((perm-x)**2)) for x in permmedoids]) + 1)
                            #dist = math.log( min( [abs(perm - x) for x in permmedoids] ) + 1)
                            distances.append(dist)
                        clump = float(sum(distances)) / float(len(permsample))
                        clumppermutations[key].append(clump)
                    except KeyError:
                        # print("keyerror", ii)
                        pass
        return clumppermutations

    def performpermutations(self,casefile,controlfile,numberofperm,afcutoff):
        '''Perform permutation testing to get a pvalue'''

        controldb,controldomains=self.fill(controlfile,afcutoff)
        #print(controldb)
        controlname=list(controldb.keys())[0][2]
        controlvalues=self.controlclump(controldb)
        controldb=self.permparse(controlfile,afcutoff)
        #print(controldb)
        casedb=self.permparse(casefile,afcutoff)
        controlpermutations=self.permutedb(controldb,numberofperm)
        casepermutations=self.permutedb(casedb,numberofperm)
        
        for record in self.clumps:
            
            protID=record[1]

            try:####added when controlpermutations[protID] gives keyerror
                controls=np.array(controlpermutations[protID])
                cases=np.array(casepermutations[protID])
                controlrecord=(record[0],record[1],controlname)
                overallresults=controlvalues[controlrecord]-self.clumps[record]
                #print(controlpermutations)
                #print(casepermutations)
                diffdist=controls-cases
                pgreater=float(len(diffdist[diffdist<overallresults]))/float(len(diffdist))
                pless=float(len(diffdist[diffdist>overallresults]))/float(len(diffdist))
                print('\t'.join(record)+'\t'+str(overallresults)+'\t'+str(pgreater)+'\t'+str(pless)+'\t'+str(self.clumps[record])+'\t'+str(controlvalues[controlrecord]))
            except KeyError:
                
                pass
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def main():

    #print("START")
    opts=CmdOpts()
    MyParser = Parser(opts.filename,opts.proteinlengths,float(opts.affreq),opts.nmutcutoff,opts.normalize, opts.structurepath)
    #MyParser()
    if(opts.controlperms and opts.permutations):
        if opts.titles:
            print('GENE\tPROTEIN_ID\tSTUDY_NAME\tCLUMP_SCORE_DIFFERENCE(Controls-Cases)\tP-value(Probability Cases have a CLUMP score greater than the Controls)\tP-value(Probability Cases have a CLUMP score less than Controls)\tCASES_Raw_Clump_Score\tCONTROL_Raw_Clump_Score')
        MyParser.performpermutations(opts.filename,opts.controlperms,int(opts.permutations),float(opts.affreq))
    else:
        if opts.titles:
            print("GENE\tPROTEIN_ID\tSTUDY_NAME\tRaw_Clump_Score")
        MyParser()
if __name__ == '__main__':
    main()

