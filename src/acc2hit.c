
/***********************************************************************
 * fComb.c
 *
 * Identify multi-hit combinations based on LR-
 *
 * Calling Parameters: Tumor matrix file
 *                     Normal gene-sample list file
 * 
 * Output: list of 3-hit combinations
 ************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/resource.h>

static const int NAME_LEN = 20;

/***********************************************************************
 *
 * get count of unique genes and samples from header of input file
 *
 ************************************************************************/
// reads first line of file, BRCA.maf2dat.matrix.out.training
void getNumGenesSamplesTumor( FILE *fp_gene_sample_matrix, int *num_genes, int *num_samples )
{
   int     i, j, ret_value;
   char    *line = NULL;
   size_t  len = 0;
   ssize_t read;

   /* First line contains number of genes and samples */
   // getline(char **line, size_t *n, FILE* filePointer)
   // Gets an entire line from the file pointer
   read = getline( &line, &len, fp_gene_sample_matrix );
   // sscanf puts the values in the string into variables i and j, so i and j will have
   // whatever the %d and %d respectively are from that line
   // 19411 208  -1 Gene  Test  --> It'll process 19411, skip white space, process 208, then stop
   ret_value = sscanf( line, "%d %d", &i, &j );

   // return value is how many items were succesfully parsed, need 2 here
   if (ret_value == 2 )
   {
      // store it in the original variables
      *num_genes   = i;
      *num_samples = j;
   }
   else
   {
      printf("ERROR: invalid input file header %d\n", ret_value);
      exit( 1 );
   }
}


/***********************************************************************
 * Load gene-sample matrix data from input file
 *
 * Calling Parameters: gene-sample matrix input file
 *                     number of genes
 *                     number of samples
 *                     gene_sample matrix (updated)
 *
 ************************************************************************/
void loadGeneSampleMatrixTumor( FILE *fp_gene_sample_matrix, int num_genes, 
   int num_samples, int *gene_sample_matrix, char *gene_id, int *tumor_samples_per_gene )
{
   
   int     i, j, k, n, ret_value;
   char    *line = NULL;
   char    *gene, *sample;
   size_t  len = 0;
   ssize_t read;

   // Malloc these to store one string, I think NAME_LEN is just the maximum length for buffering purposes?
   gene   = (char *)malloc( NAME_LEN * sizeof( char ) );
   sample = (char *)malloc( NAME_LEN * sizeof( char ) );

   /* initialize matrix */
   // Remember, when allocating tumor matrix, num_genes was rows and num_samples was columns
   /*
   tumor_matrix = (int  *)malloc( num_genes * num_samples * sizeof( int ) );
   gene_id = (char  *)malloc( num_genes * NAME_LEN * sizeof( char ));
   tumor_samples_per_gene = (int  *)malloc( num_genes * sizeof( int ) ); --> These were the initializations
   */
   for ( n = 0; n < num_genes; n++ )
   {
      tumor_samples_per_gene[n] = 0; // set these all to 0
      for ( j = 0; j < num_samples; j++ )
      {
         gene_sample_matrix[n * num_samples + j] = 0; // pointer arithmetic to set this all to 0
      }
   }

   // feof checks whether the end-of-file flag is seton the given FILE*
   // WRONG: only returns true AFTER trying to read past end of the file
   // Only getline moves the file pointer, so after getline reads the last line, 
   // feof will still return False, and getline will now read and return -1,

   // feof only reports that we ALREADY FAILED to read something
   while ( !feof( fp_gene_sample_matrix ) )
   {
      // get an entire line from the file pointer, fp_gene_sample_matrix, again, and len is updated too
      read = getline( &line, &len, fp_gene_sample_matrix );
      // Read the lines values into i, j, k as integers, and gene, sample as strings
      // i = gene identification
      // j = patient id
      // k = # of mutations
      ret_value = sscanf( line, "%d %d %d %s %s", &i, &j, &k, gene, sample ); 
      if ( ret_value == 5 )
      {
         // store the gene into geneid, which is a buffer
         strcpy( gene_id+(i*NAME_LEN), gene );
         // for this gene, for this patient, indicate the number of mutations
         gene_sample_matrix[i * num_samples + j] = k;
         if ( k > 0 )
         {
            // there was a mutation for this gene, increment
            tumor_samples_per_gene[i]++;
         }
      }
      else
      {
         printf("ERROR: reading data from input file %d\n", ret_value);
         exit( 1 );
      }
   }

   return;

}


/***********************************************************************
 *
 * get count of samples in normal gene-sample file
 *
 ************************************************************************/
int getNumSamplesNormal( FILE *fp_gene_sample_list )
{
   int     num_samples, last_sample;
   int     sample, ret_value;
   char    gene[NAME_LEN];
   char    *line = NULL;
   size_t  len = 0;
   ssize_t read;

   num_samples = last_sample = 0;
   // 
   /* read to end of file */
   // Again, error with the feof only returning True AFTER getline has done an invalid read
   // manifest_normal_normal.txt.geneSampleList is the file
   // File is in format Gene | Patient ID
   while ( !feof( fp_gene_sample_list ) )
   {
      read = getline( &line, &len, fp_gene_sample_list );
      ret_value = sscanf( line, "%s %d", gene, &sample );
      if (ret_value == 2 )
      {
         // Only want to count the number of patients so just skip if it's the same one
         if ( sample != last_sample )
         {
            num_samples++;
            last_sample = sample;
         }
      }
      else
      {
         printf("ERROR: invalid line in normal gene-sample list %d\n", ret_value);
         exit( 1 );
      }
   }
   // reset file pointer back to the start
   rewind( fp_gene_sample_list );

   // return number of samples we counted
   return( num_samples );
}


/***********************************************************************
 * Load gene-sample matrix data from input file
 *
 * Calling Parameters: gene-sample matrix input file
 *                     number of genes
 *                     number of samples
 *                     gene_sample matrix (updated)
 *
 ************************************************************************/
void loadGeneSampleMatrixNormal( FILE *fp_gene_sample_list, int num_genes, 
   int num_samples_normal, int *normal_matrix, char *gene_id, int *normal_samples_per_gene )
{
   int     j, n, ret_value;
   char    *line = NULL;
   char    gene[NAME_LEN];
   int     sample, new_sample_id, matrix_sample_index;
   size_t  len = 0;
   ssize_t read;
   int     *sample_id; /* to translate from sample# in file to sequential # */
   int     max_samples = 1000; /* for allocation of sample_id list */

   /* list of old to new (sequential excluding missing numbers) sample ids */
   // max_samples seems to be arbitrary, it is just bigger than the number of patients/samples
   sample_id = (int *)malloc( max_samples * sizeof( int ) );
   if ( sample_id == NULL )
   {
      printf( "ERROR: failed to allocate memory for normal gene_sample_matrix \n" );
     exit( 1 );
   }

   // initialize everything to -1 for now in this
   for ( j = 0; j < max_samples; j++ )
   {
      sample_id[j] = -1;
   }

   /* initialize matrix */
   // Same as initialization for loadTumor
   for ( n = 0; n < num_genes; n++ )
   {
      normal_samples_per_gene[n] = 0;
      for ( j = 0; j < num_samples_normal; j++ )
      {
         normal_matrix[n * num_samples_normal + j] = 0;
      }
   }

   new_sample_id = 0;
   // again, shouldn't use this paradigm for while loop
   while ( !feof( fp_gene_sample_list ) )
   {
      read = getline( &line, &len, fp_gene_sample_list );
      // Get the gene name, and patient ID, this is mutated
      ret_value = sscanf( line, "%s %d", gene, &sample );
      if ( ret_value == 2 )
      {
         // This logic is marking what order these came in sequentially
         if ( sample_id[sample] < 0 )
         {
            sample_id[sample] = new_sample_id;
            matrix_sample_index = new_sample_id;
            new_sample_id++;
         }
         else
         {
            // If already processed before, just get what place sequentially it came in already
            matrix_sample_index = sample_id[sample];
         }
         // Goal is to find the index that equals this gene in gene_id
         for ( n = 0; n < num_genes; n++ )
         {
            // Finds if the gene stored at index n is equal to the gene found from the file
            if ( strcmp( gene_id+(n*NAME_LEN), gene ) == 0 )
            {
               // in the matrix, for this gene (row), and the patient, set to 1
               normal_matrix[n * num_samples_normal + matrix_sample_index] = 1;
               // WHY IS THIS DOING THIS??? WHAT
               normal_samples_per_gene = 0;
               break;
            }
         }
      }
      else
      {
         printf("ERROR: reading data from normal gene-sample input file %d\n", ret_value);
         exit( 1 );
      }
   }

   free( sample_id );
}


/***********************************************************************
 * Find combination with smallest lr- 
 *
 * Calling Parameters: tumor gene-sample matrix
 *                     number of genes
 *                     number of samples
 *                     normal gene-sample matrix
 *                     tumor genes_per_sample for lr bound
 * 
 ************************************************************************/
float maxF( int *tumor_matrix, int num_genes, int num_samples_tumor, 
   int *normal_matrix, int num_samples_normal, int *tumor_samples_per_gene, 
   int *normal_samples_per_gene, int *excluded_samples, 
   int *gene1, int *gene2, float beta )
{
   int   i1, i2, j;
   int   true_pos, false_pos, true_neg, false_neg;
   float ratio, total_ratio, adj_tumor_count, adj_normal_count;
   float f;          /* f-measure */
   float f_max, f_bound, prec, recall, prec_bound, recall_bound, beta2;
   int   temp_tp, temp_fp, num_skipped;  

   beta2 = beta * beta;

   f_max = 0.0; // best score thus far
   num_skipped = 0;
   // Go through each gene combination
   for ( i1 = 0; i1 < num_genes; i1++ )
   {
      for ( i2 = i1+1; i2 < num_genes; i2++ )
      {
         // temp_tp is the minimum of the 2 genes tumor counts
         // upper bound on the maximum number of genes that have BOTH of these mutated
         temp_tp = tumor_samples_per_gene[i1];
//         temp_fp = normal_samples_per_gene[i1];
         if ( tumor_samples_per_gene[i2] < temp_tp )
         {
            temp_tp = tumor_samples_per_gene[i2];
         }
         
//         if ( normal_samples_per_gene[i2] < temp_fp )
//         {
//            temp_fp = normal_samples_per_gene[i2];
//         }
//         prec_bound   = (float) (temp_tp) / (float) (temp_tp + temp_fp);
//         recall_bound = (float) (temp_tp) / (float) num_samples_tumor;
//         f_bound = 2.0 * prec_bound * recall_bound / (prec_bound + recall_bound);
//         f_bound = (1.0 + beta2) * prec_bound * recall_bound / (beta2 * prec_bound + recall_bound);

         // Calculate the max score we can get with temp_tp
         // temp_tp is upper bound on true_pos, num_samples_normal is upper bound on true_neg
         f_bound = (float)(beta * (float) temp_tp + (float) num_samples_normal) / (float)(num_samples_tumor + num_samples_normal);

         // only process if it can beat the best we've seen before
         if ( f_bound > f_max )
         {
            true_pos = 0;
            for ( j = 0; j < num_samples_tumor; j++ )
            {
               // if this sample isn't excluded AND both genes i1 and i2 are mutated, increment true_pos
               if ( excluded_samples[j] == 0 )
               {
                  if ( ( tumor_matrix[i1 * num_samples_tumor + j] > 0 ) &&
                       ( tumor_matrix[i2 * num_samples_tumor + j] > 0 ) )
                  {
                     true_pos++;
                  }
               }
            }
            false_pos = 0;
            // false positive counts the number of times both genes are mutated in normal samples
            for ( j = 0; j < num_samples_normal; j++ )
            {
               if ( ( normal_matrix[i1 * num_samples_normal + j] > 0 ) &&
                    ( normal_matrix[i2 * num_samples_normal + j] > 0 ) )
               {
                  false_pos++;
               }
            }
            if ( true_pos > 0 ) /* avoid divide by zero */
            {
               false_neg = num_samples_tumor  - true_pos;
               true_neg  = num_samples_normal - false_pos; 
//               prec      = (float) true_pos / (float) (true_pos + false_pos);
//               recall    = (float) true_pos / (float) num_samples_tumor;
//               f         = 2.0 * prec * recall / (prec + recall );
//               f         = (1.0 + beta2) * prec * recall / (beta2 * prec + recall );
               // number of times both genes mutated in tumor + (normal samples - both mutated in normal) / (all genes together)
               f         = (float)(beta * (float) true_pos + (float) true_neg) / (float)(num_samples_tumor + num_samples_normal);
               if ( f > f_max )
               {
                  f_max  = f;
                  // set these genes for the listCombs function to use, it's the current best
                  *gene1 = i1;
                  *gene2 = i2;
               }
            }
//            else
//            {
//               printf( "True neg = %d - %d = 0", num_samples_normal, false_pos );
//            }
         }
         // otherwise just skip it
         else
         {
            num_skipped++;
         }
      }
   }

//   printf( "num skipped = %d \n", num_skipped );

   // return the best f_max score acheived
   return( f_max );
}

/***********************************************************************
 * exclude samples contatining gene combinations found
 *
 * Calling Parameters: gene index
 *                     excluded samples list (for update)
 *                     tumor gene-sample matrix
 *                     number of samples
 * 
 ************************************************************************/
int excludeSamples( int gene1, int gene2,  int *excluded_samples, int *tumor_matrix, int num_samples_tumor ) 
{
   int   s, num_excluded;

   // Exclude samples if they contain BOTH gene1 and gene2
   num_excluded = 0;
   for ( s = 0; s < num_samples_tumor; s++ )
   {
      if ( excluded_samples[s] == 0 )
      {
         if ( ( tumor_matrix[gene1 * num_samples_tumor + s] > 0 ) &&
              ( tumor_matrix[gene2 * num_samples_tumor + s] > 0 ) )
         {
            excluded_samples[s] = 1;
            num_excluded++;
         }
      }
   }
   return( num_excluded );
}

/***********************************************************************
 * list all combinations found in gene-sample matrix
 *
 * Calling Parameters: tumor gene-sample matrix 
 *                     number of genes
 *                     number of samples (tumor)
 *                     normal gene-sample matrix
 *                     number of samples (normal)
 *                     gene_id (name list)
 *                     tumor samples per gene 
 * 
 ************************************************************************/
int listCombs( int *tumor_matrix, int num_genes, int num_samples_tumor, 
   int *normal_matrix, int num_samples_normal, char *gene_id, 
   int *tumor_samples_per_gene, int *normal_samples_per_gene, float beta )
{
      /* Check combinations of genes for coverage of samples */
   // tumor_matrix[num_gene][num_patient] = number of mutations for this combination
   // num_genes = total number of genes there are
   // num_samples_tumor = total number of samples/patients there are
   // normal_matrix = normal_matrix[num_gene][num_samples_normal] = 1 if found in file, otherwise 0
   // num_samples_normal = number of normal patients
   // gene_id array is the list of gene names
   // tumor_samples_per_gene has the number of patients with a mutation in this gene, for tumor samples
   // normal_samples_per_gene should do something similar but it's being set to 0? Need to ask about that
   // beta is 0.1
   int   i, num_found, gene1, gene2, num_excluded, tot_excluded;
   float f_max;
   char *gene1_name, *gene2_name;
   int  *excluded_samples;

   gene1_name       = (char *)malloc( NAME_LEN          * sizeof( char ) );
   gene2_name       = (char *)malloc( NAME_LEN          * sizeof( char ) );
   excluded_samples = (int  *)malloc( num_samples_tumor * sizeof( int  ) );
   if ( excluded_samples == NULL )
   {
      printf( "ERROR: failed to allocate memory for excluded_samples \n" );
      exit( 1 );
   }
   // initialize this i guess
   for ( i = 0; i < num_samples_tumor; i++ )
   {
      excluded_samples[i] = 0;
   }

   num_found    = 0;
   tot_excluded = 0;
   // keep going until each has been excluded
   while ( tot_excluded < num_samples_tumor )
   {
      // Calculate the highest f_max score acheived by any 2 combinations of genes and set gene1,gene2 to be that combination
      f_max = maxF( tumor_matrix, num_genes, num_samples_tumor, normal_matrix, 
                    num_samples_normal, tumor_samples_per_gene, normal_samples_per_gene, 
                    excluded_samples, &gene1, &gene2, beta );

      // Exclude samples that include both gene1 and gene2, which affects future f_max calculations
      num_excluded = excludeSamples( gene1, gene2, excluded_samples, tumor_matrix, num_samples_tumor );
      tot_excluded += num_excluded;
      num_found++;

      // Copy the actual genen names into gene1_name and gene2_name
      strcpy( gene1_name, gene_id+(gene1*NAME_LEN) );
      strcpy( gene2_name, gene_id+(gene2*NAME_LEN) );
      printf( "%s %s %d %d F-max = %9.6f , num excluded %d, tot excluded %d \n", 
              gene1_name, gene2_name, gene1, gene2, f_max, num_excluded, tot_excluded);

      if ( num_excluded == 0 )
      {
         break;
      }
   }

   // The number of combinations found to "explain" each pair
   return( num_found );
}


/*****************************************************************************************
 * main()                                                                                *
 *                                                                                       *
 * calculate ratio of occurunce os multi-hit gene combinations                           *
 *                                                                                       *
 * Calling Parameters: list of tumor and normal gene-sample counts,                      *
 *                     list of freq mutated genes                                        *
 *                                                                                       *
 * output: ratio of occurance of multi-hit gene combinations in normal and tumor samples *
 *****************************************************************************************/
int main(int argc, char ** argv)
{
   int     num_genes, num_samples, num_samples_normal;   /* number of genes and samples in tcga maf data */
   int     num_comb;                 /* number of 3-hit combinations found in tcga samples */
   int     *tumor_matrix;            /* matrix of genesxsamples */
   int     *normal_matrix;           /* matrix of genesxsamples */
   int     *tumor_samples_per_gene;  /* total True Pos for calulating bound */
   int     *normal_samples_per_gene;  /* total False Pos for calulating bound */
   FILE    *fp_tumor_matrix;
   FILE    *fp_normal_matrix;
   char    *gene_id;                 /* list of gene ids */
   float   beta;

   // arg1 = training data, BRCA.maf2dat.matrix.out.training
   // arg2 = normal data, manifest_normal_normal.txt.training.txt.geneSampleList
   // arg3 = beta val, 0.1
   if (argc != 4)
   {
      printf("ERROR: fComb requires 2 parameters\n");
      printf("       - Tumor gene-sample count matrix file (training)\n");
      printf("       - Normal gene-sample list file (training)\n");
      printf("       - Beta value (Fbeta-score)\n");
      exit(1);
   }
   // fopen opens file, returns file pointer to read/write with
   if ( ( fp_tumor_matrix = fopen( argv[1], "r" ) ) == NULL )
   {
      printf( "ERROR: unable to open tumor gene-sample count matrix file %s, \n", argv[2] );
      exit( 1 );
   }
   // fopen the normal data too
   if ( ( fp_normal_matrix = fopen( argv[2], "r" ) ) == NULL )
   {
      printf( "ERROR: unable to open normal gene-sample count matrix file %s, \n", argv[3] );
      exit( 1 );
   }
   // atof turns string to float
   beta = atof( argv[3] );

   /* load tumor gene-sample matrix */
   // Gets the number of genes and samples, where genes was the first column and samples was second
   getNumGenesSamplesTumor( fp_tumor_matrix, &num_genes, &num_samples );
   printf( "Num Tumor genes = %d tumor samples = %d \n", num_genes, num_samples );

   // allocate memory for a matrix of [num_genes][num_samples], samples is # of patient IDs
   // Also casted the result to be an integer pointer
   tumor_matrix = (int  *)malloc( num_genes * num_samples * sizeof( int ) );
   if ( tumor_matrix == NULL )
   {
      printf( "ERROR: failed to allocate memory for tumor gene_sample_matrix \n" );
      exit( 1 );
   }

   // Name_LEN is 2, why? Casted result to be integer pointer
   gene_id = (char  *)malloc( num_genes * NAME_LEN * sizeof( char ));
   if ( gene_id == NULL )
   {
      printf( "ERROR: failed to allocate memory for gene_ids \n" );
      exit( 1 );
   }
   // Have space for an integer for each number of genes
   // For each gene, counts the number of mutations that occurs for it total across the dataset
   tumor_samples_per_gene = (int  *)malloc( num_genes * sizeof( int ) );
   if ( tumor_samples_per_gene == NULL )
   {
      printf( "ERROR: failed to allocate memory for tumor samples per gene \n" );
      exit( 1 );
   }
   // pass the allocated file pointer to the data, number of genes, number of patient samples, allocated tumor matrix, 
   // allocated list for gene ids, and allocated list for tumor samples per gene
   loadGeneSampleMatrixTumor( fp_tumor_matrix, num_genes, num_samples, 
      tumor_matrix, gene_id, tumor_samples_per_gene );

   // finished reading that into the data structures, so no longer need
   fclose( fp_tumor_matrix );

   /* load normal gene-sample matrix */
   // Just gets the number of normal samples as a simple integer
   num_samples_normal = getNumSamplesNormal( fp_normal_matrix );
   printf( "Num normal samples = %d \n", num_samples_normal );

   // Create a similar matrix as tumor_matrix with num_genes x number of patient id samples
   normal_matrix = (int *)malloc( num_genes * num_samples_normal * sizeof( int ) );
   if ( normal_matrix == NULL )
   {
      printf( "ERROR: failed to allocate memory for normal gene_sample_matrix \n" );
     exit( 1 );
   }
   
   // same thing as tumor_sample_per_gene, count the number of mutations for each gene for normal samples
   normal_samples_per_gene = (int  *)malloc( num_genes * sizeof( int ) );
   if ( normal_samples_per_gene == NULL )
   {
      printf( "ERROR: failed to allocate memory for normal samples per gene \n" );
      exit( 1 );
   }
   
   // Pass in file pointer for the normal matrix, pointing at the manifest_normal_normal.txt.geneSampleList file
   // num_genes was obtained earlier based on the tumor samples processing
   loadGeneSampleMatrixNormal( fp_normal_matrix, num_genes, num_samples_normal, 
      normal_matrix, gene_id, normal_samples_per_gene );

   // finished reading into this so done
   fclose( fp_normal_matrix );

   /* Check combinations of genes for coverage of samples */
   // MAIN LOGIC TO FIND IT HERE, EVERYTHING ELSE WAS SETUP!!
   // Reminder:
   // tumor_matrix[num_gene][num_patient] = number of mutations for this combination
   // num_genes = total number of genes there are
   // num_samples = total number of samples/patients there are
   // normal_matrix = normal_matrix[num_gene][num_samples_normal] = 1 if found in file, otherwise 0
   // number of normal patients
   // gene_id array is the list of gene names
   // tumor_samples_per_gene has the number of patients with a mutation in this gene, for tumor samples
   // normal_samples_per_gene should do something similar but it's being set to 0? Need to ask about that
   // beta is 0.1
   num_comb = listCombs( tumor_matrix, num_genes, num_samples, 
      normal_matrix, num_samples_normal, gene_id, 
      tumor_samples_per_gene, normal_samples_per_gene, beta );
   
   // Print the number of 2 hit combinations at the end
   printf( "Num 2-hit combinations = %d  (beta = %f )\n", num_comb, beta );

   // Free all memory
   free( tumor_matrix );
   free( normal_matrix );
   free( gene_id );
   free( tumor_samples_per_gene );
   free( normal_samples_per_gene );

   return( 0 );
}

