
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
 #include <omp.h>
 
 /* CSR-style index list mapping each gene to the samples where it appears */
 typedef struct {
    int *indices;   /* flattened sample indices for all genes */
    int *offsets;   /* offsets per gene into indices (length num_genes + 1) */
 } GeneSampleIndexList;
 
 static void buildSampleIndexList( const int *matrix, int num_genes, int num_samples,
    const int *samples_per_gene, GeneSampleIndexList *list );
 static void freeSampleIndexList( GeneSampleIndexList *list );
 
 
 
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
 
    #pragma omp parallel for private(j)
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
    /* Robust read: stop at EOF and tolerate trailing blank lines */
    // get an entire line from the file pointer, fp_gene_sample_matrix, again, and len is updated too
    while ( (read = getline( &line, &len, fp_gene_sample_matrix )) != -1 )
    {
       
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
    /* Robust read: stop at EOF and tolerate trailing blank lines */
    while ( (read = getline( &line, &len, fp_gene_sample_list )) != -1 )
    {
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
 
 
    #pragma omp parallel for private(j)
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
    /* Robust read: stop at EOF and tolerate trailing blank lines */
    while ( (read = getline( &line, &len, fp_gene_sample_list )) != -1 )
    {
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
                normal_samples_per_gene[n]++;
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
 
 /* Compress a dense gene x sample matrix into CSR-like sample lists per gene */
 static void buildSampleIndexList( const int *matrix, int num_genes, int num_samples,
    const int *samples_per_gene, GeneSampleIndexList *list )
 {
    int total_entries;
 
    list->indices = NULL;
    list->offsets = NULL;
 
    total_entries = 0;
    list->offsets = (int *)malloc( (num_genes + 1) * sizeof( int ) );
    if ( list->offsets == NULL )
    {
       printf( "ERROR: failed to allocate memory for sample index offsets\n" );
       exit( 1 );
    }
 
    list->offsets[0] = 0;
    for ( int g = 0; g < num_genes; g++ )
    {
       total_entries += samples_per_gene[g];
       list->offsets[g + 1] = total_entries;
    }
 
    if ( total_entries == 0 )
    {
       total_entries = 1; /* still allocate to allow pointer arithmetic */
    }
 
    list->indices = (int *)malloc( total_entries * sizeof( int ) );
    if ( list->indices == NULL )
    {
       printf( "ERROR: failed to allocate memory for sample index list\n" );
       exit( 1 );
    }
 
    for ( int g = 0; g < num_genes; g++ )
    {
       int base = g * num_samples;
       int pos  = list->offsets[g];
       for ( int s = 0; s < num_samples; s++ )
       {
          if ( matrix[base + s] > 0 )
          {
             list->indices[pos++] = s;
          }
       }
       assert( pos == list->offsets[g + 1] );
    }
 }
 
 static void freeSampleIndexList( GeneSampleIndexList *list )
 {
    free( list->indices );
    free( list->offsets );
    list->indices = NULL;
    list->offsets = NULL;
 }
 
 /* Merge-count intersection while optionally skipping excluded samples */
 static inline int countOverlapFiltered( const int *idx1, int len1, const int *idx2, int len2,
    const int *excluded_samples )
 {
    int p1, p2, count;
 
    p1 = 0;
    p2 = 0;
    count = 0;
 
    while ( (p1 < len1) && (p2 < len2) )
    {
       int s1 = idx1[p1];
       int s2 = idx2[p2];
 
       if ( s1 == s2 )
       {
          if ( (excluded_samples == NULL) || (excluded_samples[s1] == 0) )
          {
             count++;
          }
          p1++;
          p2++;
       }
       else if ( s1 < s2 )
       {
          p1++;
       }
       else
       {
          p2++;
       }
    }
 
    return( count );
 }
 
 /* Intersect two lists and mark shared samples as excluded */
 static inline int markExcludedOverlap( const int *idx1, int len1, const int *idx2, int len2,
    int *excluded_samples )
 {
    int p1, p2, num_excluded;
 
    p1 = 0;
    p2 = 0;
    num_excluded = 0;
 
    while ( (p1 < len1) && (p2 < len2) )
    {
       int s1 = idx1[p1];
       int s2 = idx2[p2];
 
       if ( s1 == s2 )
       {
          if ( excluded_samples[s1] == 0 )
          {
             excluded_samples[s1] = 1;
             num_excluded++;
          }
          p1++;
          p2++;
       }
       else if ( s1 < s2 )
       {
          p1++;
       }
       else
       {
          p2++;
       }
    }
 
    return( num_excluded );
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
 /* Evaluate best F-score pair using sparse sample lists */
 float maxF( int num_genes, int num_samples_tumor, int num_samples_normal,
    int *tumor_samples_per_gene, const GeneSampleIndexList *tumor_index_list,
    const GeneSampleIndexList *normal_index_list, int *excluded_samples,
    int *gene1, int *gene2, float beta )
 {
    int   i1, i2;
    int   true_pos, false_pos, true_neg, false_neg;
    float f;          /* f-measure */
    float f_max, f_bound;
    int   temp_tp;
 
    f_max = 0.0f; // global best score
    int best_g1 = -1;
    int best_g2 = -1;
 
    /* Parallelize search over gene pairs */
    #pragma omp parallel private(i1,i2,true_pos,false_pos,true_neg,false_neg, \
                                 f,f_bound,temp_tp) 
    {
       float local_f_max = 0.0f;   // thread-local best score
       int   local_g1    = -1;     // thread-local best gene1
       int   local_g2    = -1;     // thread-local best gene2
 
       #pragma omp for nowait
       for ( i1 = 0; i1 < num_genes; i1++ )
       {
          for ( i2 = i1+1; i2 < num_genes; i2++ )
          {
             // temp_tp is the minimum of the 2 genes tumor counts
             // upper bound on the maximum number of genes that have BOTH of these mutated
             temp_tp = tumor_samples_per_gene[i1];
             if ( tumor_samples_per_gene[i2] < temp_tp )
             {
                temp_tp = tumor_samples_per_gene[i2];
             }
 
             // Calculate the max score we can get with temp_tp
             // temp_tp is upper bound on true_pos, num_samples_normal is upper bound on true_neg
             f_bound = (float)(beta * (float) temp_tp + (float) num_samples_normal) /
                       (float)(num_samples_tumor + num_samples_normal);
 
             // only process if it can beat this thread's best so far
             if ( f_bound > local_f_max )
             {
                /* Sparse intersection of tumor samples for genes i1 and i2 */
                const int *tumor_idx1 = tumor_index_list->indices + tumor_index_list->offsets[i1];
                const int *tumor_idx2 = tumor_index_list->indices + tumor_index_list->offsets[i2];
                int len1 = tumor_index_list->offsets[i1 + 1] - tumor_index_list->offsets[i1];
                int len2 = tumor_index_list->offsets[i2 + 1] - tumor_index_list->offsets[i2];
                true_pos = countOverlapFiltered( tumor_idx1, len1, tumor_idx2, len2, excluded_samples );
 
                /* Same idea for normal samples */
                const int *normal_idx1 = normal_index_list->indices + normal_index_list->offsets[i1];
                const int *normal_idx2 = normal_index_list->indices + normal_index_list->offsets[i2];
                int nlen1 = normal_index_list->offsets[i1 + 1] - normal_index_list->offsets[i1];
                int nlen2 = normal_index_list->offsets[i2 + 1] - normal_index_list->offsets[i2];
                false_pos = countOverlapFiltered( normal_idx1, nlen1, normal_idx2, nlen2, NULL );
 
                if ( true_pos > 0 ) /* avoid divide by zero */
                {
                   false_neg = num_samples_tumor  - true_pos;
                   true_neg  = num_samples_normal - false_pos;
 
                   // number of times both genes mutated in tumor + (normal samples - both mutated in normal) / (all genes together)
                   f = (float)(beta * (float) true_pos + (float) true_neg) /
                       (float)(num_samples_tumor + num_samples_normal);
 
                   if ( f > local_f_max )
                   {
                      local_f_max = f;
                      local_g1    = i1;
                      local_g2    = i2;
                   }
                }
             }
          }
       }
 
       // Reduce per-thread bests into global best
       #pragma omp critical
       {
          if ( local_f_max > f_max )
          {
             f_max   = local_f_max;
             best_g1 = local_g1;
             best_g2 = local_g2;
          }
       }
    }
 
    *gene1 = best_g1;
    *gene2 = best_g2;
 
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
 
 
 /* Remove samples that contain both genes using their sparse lists */
 int excludeSamples( int gene1, int gene2,  int *excluded_samples, const GeneSampleIndexList *tumor_index_list ) 
 {
    const int *idx1 = tumor_index_list->indices + tumor_index_list->offsets[gene1];
    const int *idx2 = tumor_index_list->indices + tumor_index_list->offsets[gene2];
    int len1 = tumor_index_list->offsets[gene1 + 1] - tumor_index_list->offsets[gene1];
    int len2 = tumor_index_list->offsets[gene2 + 1] - tumor_index_list->offsets[gene2];
 
    return markExcludedOverlap( idx1, len1, idx2, len2, excluded_samples );
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
 int listCombs( int num_genes, int num_samples_tumor, int num_samples_normal, char *gene_id,
    int *tumor_samples_per_gene, const GeneSampleIndexList *tumor_index_list,
    const GeneSampleIndexList *normal_index_list, float beta )
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
       f_max = maxF( num_genes, num_samples_tumor, num_samples_normal,
                     tumor_samples_per_gene, tumor_index_list, normal_index_list,
                     excluded_samples, &gene1, &gene2, beta );
 
       // Exclude samples that include both gene1 and gene2, which affects future f_max calculations
       num_excluded = excludeSamples( gene1, gene2, excluded_samples, tumor_index_list );
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
    omp_set_num_threads(16);
    int     num_genes, num_samples, num_samples_normal;   /* number of genes and samples in tcga maf data */
    int     num_comb;                 /* number of 3-hit combinations found in tcga samples */
    int     *tumor_matrix;            /* matrix of genesxsamples */
    int     *normal_matrix;           /* matrix of genesxsamples */
    int     *tumor_samples_per_gene;  /* total True Pos for calulating bound */
    int     *normal_samples_per_gene;  /* total False Pos for calulating bound */
    GeneSampleIndexList tumor_index_list;
    GeneSampleIndexList normal_index_list;
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
       printf( "ERROR: unable to open tumor gene-sample count matrix file %s, \n", argv[1] );
       exit( 1 );
    }
    // fopen the normal data too
    if ( ( fp_normal_matrix = fopen( argv[2], "r" ) ) == NULL )
    {
       printf( "ERROR: unable to open normal gene-sample count matrix file %s, \n", argv[2] );
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
 
    buildSampleIndexList( tumor_matrix, num_genes, num_samples, 
       tumor_samples_per_gene, &tumor_index_list );
    buildSampleIndexList( normal_matrix, num_genes, num_samples_normal, 
       normal_samples_per_gene, &normal_index_list );
 
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
    num_comb = listCombs( num_genes, num_samples, num_samples_normal, gene_id, 
       tumor_samples_per_gene, &tumor_index_list, &normal_index_list, beta );
    
    // Print the number of 2 hit combinations at the end
    printf( "Num 2-hit combinations = %d  (beta = %f )\n", num_comb, beta );
 
    // Free all memory
    freeSampleIndexList( &tumor_index_list );
    freeSampleIndexList( &normal_index_list );
    free( tumor_matrix );
    free( normal_matrix );
    free( gene_id );
    free( tumor_samples_per_gene );
    free( normal_samples_per_gene );
 
    return( 0 );
 }