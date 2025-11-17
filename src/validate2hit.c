
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
         gene_sample_matrix[i * num_samples + j] = k;
         gene[NAME_LEN - 1] = '\0';
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
   while ( !feof( fp_comb ) )
   {
      read = getline( &line, &len, fp_comb );
      n++;
   }

   return( n - 1 ); /* one extra read before eof detected */

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
   for ( j = 0; j < num_samples; j++ )
   {
      num_comb_per_sample = 0;
      for ( i = 0; i < num_comb; i++ )
      {
         g1 = comb_list[i * 2 + 0];
         g2 = comb_list[i * 2 + 1];

         if ( ( g1 > -1 ) && ( g2 > -1 ) )
         {
            if ( ( gene_sample_matrix[g1 * num_samples + j] > 0 ) &&
                 ( gene_sample_matrix[g2 * num_samples + j] > 0 ) )
            {
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

   if ( ( fp_gene_sample_matrix = fopen( argv[1], "r" ) ) == NULL )
   {
      printf( "ERROR: geneComb cannot open gene-sample count matrix file %s, \n", argv[1] );
      exit( 1 );
   }

   if ( ( fp_comb = fopen( argv[2], "r" ) ) == NULL )
   {
      printf( "ERROR: geneComb cannot open combination list file %s, \n", argv[1] );
      exit( 1 );
   }

   /* read maf2dat file to get number of gene and sample ids */
   getNumGenesSamples( fp_gene_sample_matrix, &num_genes, &num_samples );
   printf("Num genes = %d samples = %d \n", num_genes, num_samples);

   /* allocate space and load matrix data */
   gene_sample_matrix = (int  *)malloc( num_genes   * num_samples * sizeof( int ));
   gene_list          = (char *)malloc( num_genes   * 20 * sizeof( char ));
   loadGeneSampleMatrix( fp_gene_sample_matrix, num_genes, num_samples, 
      gene_sample_matrix, gene_list );
   fclose( fp_gene_sample_matrix );

   /* get number of combinations */
   num_comb = getNumComb( fp_comb );
   printf("Num combinations = %d\n", num_comb);

   /* alloate space and load list of combinations */
   comb_list = (int *) malloc( num_comb * 2 * sizeof( int ) );
   rewind( fp_comb );
   loadComb( num_genes, gene_list, num_comb, fp_comb, comb_list );
   fclose( fp_comb );

   percent_found = countCombPerSample( num_genes, num_samples, gene_sample_matrix, 
      num_comb, comb_list );
   printf( "Percent found %f \n", percent_found );

   free( gene_sample_matrix );
   free( gene_list );
//   free( comb_list );

   return( 0 );
}

