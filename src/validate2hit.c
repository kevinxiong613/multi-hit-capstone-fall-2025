
/***********************************************************************
 * geneComb.c
 *
 * Determine best set of 3-hit combinations
 *
 * Calling Parameters: gene-sample mutation matrix, list of possible combindations
 * 
 * return: None.
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
void getNumGenesSamples( FILE *fp_gene_sample_matrix, int *num_genes, int *num_samples )
{
   int     i, j, ret_value;
   char    *line = NULL;
   size_t  len = 0;
   ssize_t read;

   /* First line contains number of genes and samples */
   read = getline( &line, &len, fp_gene_sample_matrix );
   // Process the top 2 lines and get the number of genes and number of samples
   ret_value = sscanf( line, "%d %d", &i, &j );
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
 * Calling Parameters: Gene-sample matrix input file
 *                     Updateable matrix
 * 
 ************************************************************************/
void loadGeneSampleMatrix( FILE *fp_gene_sample_matrix, int num_genes, 
   int num_samples, int *gene_sample_matrix, char *gene_list )
{
   int     i, j, k, ret_value;
   char    *line = NULL;
   char    gene[NAME_LEN], sample[NAME_LEN];
   size_t  len = 0;
   ssize_t read;


   while ( !feof( fp_gene_sample_matrix ) )
   {
      read = getline( &line, &len, fp_gene_sample_matrix );
      ret_value = sscanf( line, "%d %d %d %s %s", &i, &j, &k, gene, sample );
      if (ret_value == 5 )
      {
         // Set the number of mutations this tumor matrix has
         gene_sample_matrix[i * num_samples + j] = k;
         // For some reason here we set a terminating character too
         gene[NAME_LEN - 1] = '\0';
         // Copy this into the gene list based on the ID
         strcpy( gene_list + (i * NAME_LEN), gene );
      }
      else
      {
         printf("ERROR: reading data from input file %d\n", ret_value);
         exit( 1 );
      }
   }
}

/***********************************************************************
 * get number of matched gene ids between TCGA and COSMIC
 *
 * Calling Parameters: number of genes in TCGA lsit
*                      TCGA genes list
 *                     COSMIC file pointer
 *                     number of matched COSMIC genex 
 * 
 ************************************************************************/
int getNumComb( FILE *fp_comb )
{
   int  n;
   char *line = NULL;
   size_t  len = 0;
   ssize_t read;

   n = 0;
   // Reads through BRCA-combinations and sees how many lines there are
   while ( !feof( fp_comb ) )
   {
      read = getline( &line, &len, fp_comb );
      n++;
   }

   return( n - 1 ); /* one extra read before eof detected at*/

}

/***********************************************************************
 * load gene list from COSMIC and match to gene-list from TCGA
 *
 * Calling Parameters: TCGA genes list
 *                     COSMIC file pointer
 *                     COSMIC genes list (update)
 * 
 ************************************************************************/
void loadComb( int num_genes, char *gene_list, int num_comb, FILE *fp_comb, int *comb_list )
{
   int  i, j, g1, g2, ret_value;
   char gg1[NAME_LEN], gg2[NAME_LEN];
   char *line = NULL;
   size_t len = 0;
   ssize_t read;

   /* read data lines */
   for ( i = 0; i < num_comb; i++ )
   {
      read = getline( &line, &len, fp_comb );
      ret_value = sscanf( line, "%s %s", gg1, gg2 );
      if ( ret_value == 2 )
      {
         g1 = g2 = -1;
         // Find the indices that these two genes occur at in the gene_list
         for ( j = 0; j < num_genes; j++ )
         {
            if ( strcmp( gene_list + j * NAME_LEN, gg1 ) == 0 )
            {
               g1 = j;
            }
            if ( strcmp( gene_list + j * NAME_LEN, gg2 ) == 0 )
            {
               g2 = j;
            }
            if ( ( g1 > -1 ) && ( g2 > -1 ) )
            {
               break;
            }
         }
         // In the combination list, set the indices they occur at in the gene list
         comb_list[i * 2 + 0] = g1;
         comb_list[i * 2 + 1] = g2;
//         printf("Combination: %d %d %d %s %s %s\n", g1, g2, g3, gg1, gg2, gg3 );
      }
      else
      {
         printf("ERROR: reading data from comb list file - %d\n", ret_value);
      }
   }
}

/***********************************************************************
 * Count number of combinations for each sample
 *
 ************************************************************************/
float countCombPerSample( int num_genes, int num_samples, 
   int *gene_sample_matrix, int num_comb, int *comb_list )
{
   int g1, g2, i, j;
   int num_comb_per_sample, found, not_found;
   float percent_found;

   found = not_found = 0;
   // Go through each sample
   for ( j = 0; j < num_samples; j++ )
   {
      num_comb_per_sample = 0;
      // Go through each combination we got for acc2hit
      for ( i = 0; i < num_comb; i++ )
      {
         // We can get the indices they occur at since we saved that in the comb_list
         g1 = comb_list[i * 2 + 0];
         g2 = comb_list[i * 2 + 1];

         if ( ( g1 > -1 ) && ( g2 > -1 ) )
         {
            // If both are mutated for this sample
            if ( ( gene_sample_matrix[g1 * num_samples + j] > 0 ) &&
                 ( gene_sample_matrix[g2 * num_samples + j] > 0 ) )
            {
               // Increment this
               num_comb_per_sample++;
//               printf("found combination %d in sample %d\n", i, j );
            }
         }
      }
      printf("%5d combs for sample %5d\n", num_comb_per_sample, j);
      if ( num_comb_per_sample > 0 )
      {
         found++;
      }
      else
      {
         not_found++;
      }
   }
   percent_found = (float) found / (float) ( found + not_found );
   return( percent_found );
}

/***********************************************************************
 * main()
 *
 * identify 3-hit coimbinations
 *
 ************************************************************************/
int main(int argc, char ** argv)
{
   int     num_genes, num_samples;  /* number of genes and samples in tcga maf data */
   int     num_comb;                /* number of 3-hit combinations found in tcga samples */
   int     *gene_sample_matrix;     /* matrix of genesxsamples */
   char    *gene_list;              /* list of gene names */
   int     *comb_list;              /* list of matching gene-index in cosmic */
   FILE    *fp_gene_sample_matrix, *fp_comb;
   float   percent_found;

   /* validate and open files */
   if (argc != 3)
   {
      printf("ERROR: geneComb requires 2 parameters - gene-sample count matrix file and combination list\n");
      exit(1);
   }
   // BRCA.maf2dat.matrix.out.test is the file
   if ( ( fp_gene_sample_matrix = fopen( argv[1], "r" ) ) == NULL )
   {
      printf( "ERROR: geneComb cannot open gene-sample count matrix file %s, \n", argv[1] );
      exit( 1 );
   }

   // result/$cancer/$cancer-combinations file
   if ( ( fp_comb = fopen( argv[2], "r" ) ) == NULL )
   {
      printf( "ERROR: geneComb cannot open combination list file %s, \n", argv[1] );
      exit( 1 );
   }

   /* read maf2dat file to get number of gene and sample ids */
   // Again, just reads the top 2 lines
   getNumGenesSamples( fp_gene_sample_matrix, &num_genes, &num_samples );
   printf("Num genes = %d samples = %d \n", num_genes, num_samples);

   /* allocate space and load matrix data */
   // This is the matrix for gene + patient
   gene_sample_matrix = (int  *)malloc( num_genes   * num_samples * sizeof( int ));
   // This is equivalent to gene_id from acc2hit, stores the name of genes
   gene_list          = (char *)malloc( num_genes   * 20 * sizeof( char ));
   loadGeneSampleMatrix( fp_gene_sample_matrix, num_genes, num_samples, 
      gene_sample_matrix, gene_list );
   fclose( fp_gene_sample_matrix );

   /* get number of combinations */
   // Reads through BRCA-combinations and sees how many lines there are
   num_comb = getNumComb( fp_comb );
   printf("Num combinations = %d\n", num_comb);

   /* alloate space and load list of combinations */
   comb_list = (int *) malloc( num_comb * 2 * sizeof( int ) );
   rewind( fp_comb );
   loadComb( num_genes, gene_list, num_comb, fp_comb, comb_list );
   fclose( fp_comb );

   // comb_list holds [[index of g1 in gene_list, index of g2 in gene_list], etc...]
   percent_found = countCombPerSample( num_genes, num_samples, gene_sample_matrix, 
      num_comb, comb_list );
   printf( "Percent found %f \n", percent_found );

   free( gene_sample_matrix );
   free( gene_list );
//   free( comb_list );

   return( 0 );
}

