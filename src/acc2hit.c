
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

   gene   = (char *)malloc( NAME_LEN * sizeof( char ) );
   sample = (char *)malloc( NAME_LEN * sizeof( char ) );

   /* initialize matrix */
   for ( n = 0; n < num_genes; n++ )
   {
      tumor_samples_per_gene[n] = 0;
      for ( j = 0; j < num_samples; j++ )
      {
         gene_sample_matrix[n * num_samples + j] = 0;
      }
   }

   while ( !feof( fp_gene_sample_matrix ) )
   {
      read = getline( &line, &len, fp_gene_sample_matrix );
      ret_value = sscanf( line, "%d %d %d %s %s", &i, &j, &k, gene, sample );
      if ( ret_value == 5 )
      {
         strcpy( gene_id+(i*NAME_LEN), gene );
         gene_sample_matrix[i * num_samples + j] = k;
         if ( k > 0 )
         {
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

   /* read to end of file */
   while ( !feof( fp_gene_sample_list ) )
   {
      read = getline( &line, &len, fp_gene_sample_list );
      ret_value = sscanf( line, "%s %d", gene, &sample );
      if (ret_value == 2 )
      {
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
   rewind( fp_gene_sample_list );

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
   sample_id = (int *)malloc( max_samples * sizeof( int ) );
   if ( sample_id == NULL )
   {
      printf( "ERROR: failed to allocate memory for normal gene_sample_matrix \n" );
     exit( 1 );
   }
   for ( j = 0; j < max_samples; j++ )
   {
      sample_id[j] = -1;
   }

   /* initialize matrix */
   for ( n = 0; n < num_genes; n++ )
   {
      normal_samples_per_gene[n] = 0;
      for ( j = 0; j < num_samples_normal; j++ )
      {
         normal_matrix[n * num_samples_normal + j] = 0;
      }
   }

   new_sample_id = 0;
   while ( !feof( fp_gene_sample_list ) )
   {
      read = getline( &line, &len, fp_gene_sample_list );
      ret_value = sscanf( line, "%s %d", gene, &sample );
      if ( ret_value == 2 )
      {
         if ( sample_id[sample] < 0 )
         {
            sample_id[sample] = new_sample_id;
            matrix_sample_index = new_sample_id;
            new_sample_id++;
         }
         else
         {
            matrix_sample_index = sample_id[sample];
         }
         for ( n = 0; n < num_genes; n++ )
         {
            if ( strcmp( gene_id+(n*NAME_LEN), gene ) == 0 )
            {
               normal_matrix[n * num_samples_normal + matrix_sample_index] = 1;
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

   f_max = 0.0;
   num_skipped = 0;
   for ( i1 = 0; i1 < num_genes; i1++ )
   {
      for ( i2 = i1+1; i2 < num_genes; i2++ )
      {
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
         f_bound = (float)(beta * (float) temp_tp + (float) num_samples_normal) / (float)(num_samples_tumor + num_samples_normal);
         if ( f_bound > f_max )
         {
            true_pos = 0;
            for ( j = 0; j < num_samples_tumor; j++ )
            {
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
               f         = (float)(beta * (float) true_pos + (float) true_neg) / (float)(num_samples_tumor + num_samples_normal);
               if ( f > f_max )
               {
                  f_max  = f;
                  *gene1 = i1;
                  *gene2 = i2;
               }
            }
//            else
//            {
//               printf( "True neg = %d - %d = 0", num_samples_normal, false_pos );
//            }
         }
         else
         {
            num_skipped++;
         }
      }
   }

//   printf( "num skipped = %d \n", num_skipped );

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
   for ( i = 0; i < num_samples_tumor; i++ )
   {
      excluded_samples[i] = 0;
   }

   num_found    = 0;
   tot_excluded = 0;
   while ( tot_excluded < num_samples_tumor )
   {
      f_max = maxF( tumor_matrix, num_genes, num_samples_tumor, normal_matrix, 
                    num_samples_normal, tumor_samples_per_gene, normal_samples_per_gene, 
                    excluded_samples, &gene1, &gene2, beta );

      num_excluded = excludeSamples( gene1, gene2, excluded_samples, tumor_matrix, num_samples_tumor );
      tot_excluded += num_excluded;
      num_found++;

      strcpy( gene1_name, gene_id+(gene1*NAME_LEN) );
      strcpy( gene2_name, gene_id+(gene2*NAME_LEN) );
      printf( "%s %s %d %d F-max = %9.6f , num excluded %d, tot excluded %d \n", 
              gene1_name, gene2_name, gene1, gene2, f_max, num_excluded, tot_excluded);

      if ( num_excluded == 0 )
      {
         break;
      }
   }


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
   getNumGenesSamplesTumor( fp_tumor_matrix, &num_genes, &num_samples );
   printf( "Num Tumor genes = %d tumor samples = %d \n", num_genes, num_samples );

   tumor_matrix = (int  *)malloc( num_genes * num_samples * sizeof( int ) );
   if ( tumor_matrix == NULL )
   {
      printf( "ERROR: failed to allocate memory for tumor gene_sample_matrix \n" );
      exit( 1 );
   }
   gene_id = (char  *)malloc( num_genes * NAME_LEN * sizeof( char ));
   if ( gene_id == NULL )
   {
      printf( "ERROR: failed to allocate memory for gene_ids \n" );
      exit( 1 );
   }
   tumor_samples_per_gene = (int  *)malloc( num_genes * sizeof( int ) );
   if ( tumor_samples_per_gene == NULL )
   {
      printf( "ERROR: failed to allocate memory for tumor samples per gene \n" );
      exit( 1 );
   }
   loadGeneSampleMatrixTumor( fp_tumor_matrix, num_genes, num_samples, 
      tumor_matrix, gene_id, tumor_samples_per_gene );

   fclose( fp_tumor_matrix );

   /* load normal gene-sample matrix */
   num_samples_normal = getNumSamplesNormal( fp_normal_matrix );
   printf( "Num normal samples = %d \n", num_samples_normal );

   normal_matrix = (int *)malloc( num_genes * num_samples_normal * sizeof( int ) );
   if ( normal_matrix == NULL )
   {
      printf( "ERROR: failed to allocate memory for normal gene_sample_matrix \n" );
     exit( 1 );
   }
   normal_samples_per_gene = (int  *)malloc( num_genes * sizeof( int ) );
   if ( normal_samples_per_gene == NULL )
   {
      printf( "ERROR: failed to allocate memory for normal samples per gene \n" );
      exit( 1 );
   }
   loadGeneSampleMatrixNormal( fp_normal_matrix, num_genes, num_samples_normal, 
      normal_matrix, gene_id, normal_samples_per_gene );

   fclose( fp_normal_matrix );

/* Check combinations of genes for coverage of samples */
   num_comb = listCombs( tumor_matrix, num_genes, num_samples, 
      normal_matrix, num_samples_normal, gene_id, 
      tumor_samples_per_gene, normal_samples_per_gene, beta );
   printf( "Num 2-hit combinations = %d  (beta = %f )\n", num_comb, beta );

   free( tumor_matrix );
   free( normal_matrix );
   free( gene_id );
   free( tumor_samples_per_gene );
   free( normal_samples_per_gene );

   return( 0 );
}

