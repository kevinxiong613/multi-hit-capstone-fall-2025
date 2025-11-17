#! /usr/bin/python


# get list of mut gene mut in sample
def countGeneMut(samp_file, comb_file):

   genes       = []
   num_genes   = 0

#  get number of mut gene
   with open(samp_file, 'r') as samp:

      for line in samp:
         fields = line.split() # space separated

         gene_id = fields[0]
         if (gene_id not in genes):
            genes.append(gene_id)
            num_genes += 1

   samp.close()
   print 'Number of mutated genes:', num_genes

   with open(comb_file, 'r') as comb:

      num_found = 0
      for line in comb:
         fields = line.split() # space separated

         gene1 = fields[0]
         gene2 = fields[1]
         if ( (gene1 in genes) and (gene2 in genes) ):
            print gene1, gene2
            num_found += 1

   comb.close()
   print 'Number of combinations found:', num_found


   return 

#
# main routine for validateNormal.py
#

if __name__ == "__main__":
    import sys
    import math
    import numpy as np


    if (len(sys.argv) < 3):
       print 'validateNormal.py requires 2 paramenters:'
       print ' - Normal sample mut gene list file'
       print ' - Combinations file'
       exit(1) 

    samp_file  = sys.argv[1] 
    comb_file  = sys.argv[2] 
    print 'Validating Normal sample:', samp_file
    print '   using:', comb_file

    countGeneMut(samp_file, comb_file)

    exit(0)

