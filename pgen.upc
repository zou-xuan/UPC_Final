#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <upc_io.h>

#include "packingDNAseq.h"
#include "contig_generation.h"

typedef struct s_kmer_t s_kmer_t;
struct s_kmer_t{
   char kmer[KMER_PACKED_LENGTH];
   char l_ext;
   char r_ext;
   int64_t next; // Store the index of next, to save us from the messy pointer stuff
};

/* Auxiliary function for computing hash values */
int64_t hashseq(int64_t  hashtable_size, char *seq, int size)
{
   unsigned long hashval;
   hashval = 5381;
   for(int i = 0; i < size; i++) {
      hashval = seq[i] +  (hashval << 5) + hashval;
   }

   return hashval % hashtable_size;
}

/* Returns the hash value of a kmer */
int64_t hashkmer(int64_t  hashtable_size, char *seq)
{
   return hashseq(hashtable_size, seq, KMER_PACKED_LENGTH);
}

// the edited add_kmer function 
void s_add_kmer( shared int64_t* hashtable, int tablesize, 
					shared s_kmer_t * memory_heap, int64_t heap_pos,
					const unsigned char *kmer, char left_ext, char right_ext)
{
	/* Pack a k-mer sequence appropriately */
	char packedKmer[KMER_PACKED_LENGTH];
	packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
	int64_t hashval = hashkmer(tablesize, (char*) packedKmer);

	s_kmer_t tmp;

	memcpy(tmp.kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
	tmp.l_ext = left_ext;
	tmp.r_ext = right_ext;
	upc_memput(&memory_heap[heap_pos], &tmp, sizeof(s_kmer_t)); // Put the local variable at the shared memory

	int64_t old = -1;
	int64_t next_one = old;
	do{
		memory_heap[heap_pos].next = next_one;
		old = next_one;
		next_one = bupc_atomicI64_cswap_relaxed(&hashtable[hashval], old, heap_pos);
	} while (next_one != old);
	// Effectively adding to the begining of the linked list of the specified bucket
	// while loop and atomic operation prevents the race condition and makes sure the current index is inserted
}


/* Looks up a kmer in the hash table and returns a pointer to that entry */
int64_t s_lookup_kmer(shared int64_t* hashtable, int tablesize, shared s_kmer_t * memory_heap, const unsigned char *kmer)
{
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(tablesize, (char*) packedKmer);

   int kmer_idx = hashtable[hashval];
   s_kmer_t tmp = memory_heap[kmer_idx];

   for (; kmer_idx!=-1; ) {
      if ( memcmp(packedKmer, tmp.kmer, KMER_PACKED_LENGTH * sizeof(char)) == 0 ) {
         return kmer_idx;
      }
      kmer_idx = tmp.next;
      tmp = memory_heap[kmer_idx];
   }
   return -1; // Should not reach here
}


int main(int argc, char *argv[]){

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;

	char cur_contig[MAXIMUM_CONTIG_SIZE], unpackedKmer[KMER_LENGTH+1], left_ext, right_ext;
	int64_t posInContig, contigID = 0, totBases = 0;
	unpackedKmer[KMER_LENGTH] = '\0';

	/** Read input **/
	upc_barrier;
	inputTime -= gettime();

	///////////////////////////////////////////
	// Your code for input file reading here //
	///////////////////////////////////////////
	/* Initialize lookup table that will be used for the DNA packing routines */
	init_LookupTable();

	/* Extract the number of k-mers in the input file */
	/* Read the input file name */
	char* input_UFX_name = argv[1];
	int64_t nKmers = getNumKmersInUFX(input_UFX_name);

	/* Compute the nKmers for each thread */
	int64_t nKmers_per_thread = (nKmers + THREADS - 1) / THREADS; // The Ceiling
	int64_t nKmers_local = nKmers_per_thread;
	if(MYTHREAD*nKmers_per_thread >= nKmers) nKmers_local = 0; // The last several threads are allocated 0 kmer
	else if((MYTHREAD+1)*nKmers_per_thread >= nKmers) nKmers_local = nKmers - MYTHREAD*nKmers_per_thread;

	int64_t nKmers_before = MYTHREAD * nKmers_per_thread;
	if(nKmers_local==0) nKmers_before = nKmers;

	int64_t chars_read = nKmers_local * LINE_SIZE;
	int64_t chars_before = nKmers_before * LINE_SIZE;

	unsigned char* working_buffer =  (unsigned char*) malloc(chars_read * sizeof(unsigned char));

	/* Read the kmers from the input file and store them in the working_buffer */
	upc_file_t *input_file;
	input_file = upc_all_fopen(input_UFX_name, UPC_RDONLY | UPC_INDIVIDUAL_FP, 0, NULL);
	upc_all_fseek(input_file, chars_before * sizeof(unsigned char), UPC_SEEK_SET);
	int64_t cur_chars_read = upc_all_fread_local(input_file, working_buffer, sizeof(unsigned char), chars_read,
		UPC_IN_ALLSYNC | UPC_OUT_ALLSYNC);

	upc_barrier;
	upc_all_fclose(input_file);

	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();


	///////////////////////////////////////////
	// Your code for graph construction here //
	///////////////////////////////////////////
	//create_hash_table(Kmer_local,&memory_heap,hashtable);

	int64_t n_buckets = nKmers * LOAD_FACTOR;
	shared int64_t *hashtable = (shared int64_t *) upc_all_alloc(n_buckets, sizeof(int64_t));

	if (hashtable == NULL) {
	  fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %lld buckets of %lu bytes\n", n_buckets, sizeof(int64_t));
	  exit(1);
	}

	int64_t i;
	upc_forall(i=0; i<n_buckets; i++; &hashtable[i])
		hashtable[i] = -1;

	int64_t heap_pos = 0; // it's like the memory_heap->posInHeap
	shared s_kmer_t* memory_heap = (shared s_kmer_t *) upc_all_alloc(THREADS, nKmers_per_thread * sizeof(s_kmer_t));

	if (memory_heap == NULL) {
	  fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
	  exit(1);
	}

	int64_t start_list_size = 0;
	int64_t * start_list = (int64_t *) malloc (nKmers_local * sizeof(int64_t));
	// Allocating nKmers_local space for the start_list, but definitely a waste!

	upc_barrier;

	int64_t ptr = 0;
	while (ptr < cur_chars_read) {
		/* working_buffer[ptr] is the start of the current k-mer                */
		/* so current left extension is at working_buffer[ptr+KMER_LENGTH+1]    */
		/* and current right extension is at working_buffer[ptr+KMER_LENGTH+2]  */

		left_ext = (char) working_buffer[ptr + KMER_LENGTH + 1];
		right_ext = (char) working_buffer[ptr + KMER_LENGTH + 2];

		/* Add k-mer to hash table */
		s_add_kmer(hashtable, n_buckets, memory_heap, heap_pos + nKmers_before,
			&working_buffer[ptr], left_ext, right_ext);

		/* Create also a list with the "start" kmers: nodes with F as left (backward) extension */
		if (left_ext == 'F') {
			start_list[start_list_size] = heap_pos + nKmers_before;
			start_list_size++;
		}

		/* Move to the next k-mer in the input working_buffer */
		ptr += LINE_SIZE;
		heap_pos++;
	}


	upc_barrier;
	constrTime += gettime();

	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////

	char filename[20];
  	sprintf (filename, "pgen%d.out", MYTHREAD);

  	FILE* output_file = fopen(filename, "w");

  	int64_t cur_kmer_idx;
	/* Pick start nodes from the startKmersList */
	for(int64_t curStartNode = 0; curStartNode < start_list_size; curStartNode++){
		/* Need to unpack the seed first */
		cur_kmer_idx = start_list[curStartNode];
		s_kmer_t cur_kmer = memory_heap[cur_kmer_idx];
		unpackSequence((unsigned char*) cur_kmer.kmer,  (unsigned char*) unpackedKmer, KMER_LENGTH);

		/* Initialize current contig with the seed content */
		memcpy(cur_contig ,unpackedKmer, KMER_LENGTH * sizeof(char));
		posInContig = KMER_LENGTH;
		right_ext = cur_kmer.r_ext;

		/* Keep adding bases while not finding a terminal node */
		while (right_ext != 'F') {
			cur_contig[posInContig] = right_ext;
			posInContig++;
			/* At position cur_contig[posInContig-KMER_LENGTH] starts the last k-mer in the current contig */
			cur_kmer_idx = s_lookup_kmer(hashtable, n_buckets, memory_heap, 
				(const unsigned char *) &cur_contig[posInContig-KMER_LENGTH]);
			cur_kmer = memory_heap[cur_kmer_idx];
			right_ext = cur_kmer.r_ext;
		}

		/* Print the contig since we have found the corresponding terminal node */
		cur_contig[posInContig] = '\0';
		fprintf(output_file,"%s\n", cur_contig);
		contigID++;
		totBases += strlen(cur_contig);
	}
	fclose(output_file);
	// printf("Generated %lld contigs with %lld total bases\n", contigID, totBases);


	upc_barrier;
	traversalTime += gettime();

	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
