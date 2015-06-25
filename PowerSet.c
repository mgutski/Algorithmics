//----------------------------------------------------------------------------
//
//	Copyright 2015 Michael Gutowski
//  	Licensed under the Apache License, Version 2.0 (the "License");
//  	you may not use this file except in compliance with the License.
//  	You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
//  	Unless required by applicable law or agreed to in writing, software
//  	distributed under the License is distributed on an "AS IS" BASIS,
//  	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  	See the License for the specific language governing permissions and
//  	limitations under the License.
//
//	The program generates a power set P(A) of the set A supplied via the 
//	standard input (stdin). A, n and k are input parameters where A is the 
//	input set, n (>=0) is the cardinality of set A and k (>=1) is the level 
//	of the power set to be generated. For example the k-th power set of set 
//	A, where k is set to 2, is P(P(A)).
//  	The resultant power set P(A) is printed to the standard output (stdout).
//	The time complexity of the algorithm used for generating the power set is 
//	O(2^n * n * k). The problem of generating the power set of a set is a  
//	decision 0-1 problem classified as NP-complete.
//	The program creates a memory pool to reduce the memory 
//	fragmentation and allow fast dynamic memory allocation.
//
//----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
//#include <malloc.h>

#if defined(WIN32) || defined(_WIN32)
__forceinline int fgetc_unlocked(FILE *stream) { return _fgetc_nolock(stream); }
#endif

//#define ONLINE_JUDGE

#define FLUSH_BUFF() while(fgetc_unlocked(stdin) != '\n');
#define OPENBRACE '{'
#define CLOSEBRACE '}'
#define CREATEEMPTYSET(_set) {											\
	_set = (char*) pool_alloc(3 * sizeof(char));						\
	*_set = OPENBRACE; *(_set+1) = CLOSEBRACE; *(_set+2) = '\0';		\
}

/***** Function Prototypes *****/
void read_args(unsigned *k, unsigned *n, char ***set);
unsigned read_uint(void);
char** read_set(const unsigned card);

unsigned powerset(char ***pset, char **set, unsigned card, unsigned k);
unsigned powerset_core(char ***pset, char **set, unsigned card);
void set_print(const char **set, const unsigned card);
void set_destroy(char **set, const unsigned card);

void pool_create(const size_t size);
void pool_create_maxalloc(void);
void pool_destroy(void);
void* pool_alloc(const size_t size);
void pool_free(void* mem);
size_t pool_rsv(void);
void* pool_rls(void);
void pool_write(char* c);

/**
 ** Program entry point.
 **/
int main(void)
{
	unsigned k, n, pcard;
	char **set = NULL, **pset = NULL;

	#ifndef ONLINE_JUDGE
	printf("[MemPool -> ALLOC]\n");
	pool_create_maxalloc();
	printf("[MemPool -> OK]\n\n");
	#else
	//pool_create(0x400); // 1KB
	pool_create(0x8000000); // 128MB
	//pool_create(0x10000000); // 256MB
	//pool_create(_HEAP_MAXREQ);
	#endif
	
	read_args(&k, &n, &set);

	#ifndef ONLINE_JUDGE
	assert(k >= 1 && n >= 0);
	#endif

	pcard = powerset(&pset, set, n, k);

	#ifndef ONLINE_JUDGE
	printf("\nP(A) = ");
	set_print(pset, pcard);
	printf("\n\n");

	system("PAUSE");
	#else 
	set_print(pset, pcard);
	#endif

	pool_destroy();

	return 0;
}

/**
 ** Reads the input parameters from the standard input stream (stdin).
 ** The data is read in order of function parameters.
 **/
void read_args(unsigned *k, unsigned *n, char ***set)
{
	#ifndef ONLINE_JUDGE
	printf("k = ");
	#endif
	
	*k = read_uint();

	#ifndef ONLINE_JUDGE
	FLUSH_BUFF();
	printf("n = ");
	#endif

	*n = read_uint();

	#ifndef ONLINE_JUDGE
	FLUSH_BUFF();
	printf("A = ");
	#endif

	*set = read_set(*n);
}

/**
 ** Reads the 32-bit unsigned integer from the standard input stream (stdin).
 **/
unsigned read_uint(void)
{
	unsigned i = 0;
	register int c = fgetc_unlocked(stdin);
	
	while(c >= '0' && c <= '9')
	{
		i = i * 10 + (c - '0');
		c = fgetc_unlocked(stdin);
	}

	return i;
}

/**
 ** Reads the set of specified cardinality.
 ** The set supplied via the standard input stream (stdin) must be a finite 
 ** subset of natural numbers.
 **/
char** read_set(const unsigned card)
{
	char c; int ci;
	unsigned i = 0;
	size_t buffsize;

	char** set = (char**) pool_alloc(card * sizeof(char*));

	while((ci = fgetc_unlocked(stdin)) != '\n' && ci != EOF && i < card)
	{
		if(ci < '0' || ci > '9') continue;
		
		buffsize = pool_rsv();
		
		for(; ci >= '0' && ci <= '9' && buffsize > 1; ci = fgetc_unlocked(stdin), buffsize--)
		{
			c = ci;
			pool_write(&c);
		}

		c = '\0';
		pool_write(&c);

		set[i++] = (char*) pool_rls();
	}
	
	return set;
}

/**
 ** Generates a k-th power set of the specified input set.
 ** The parameter 'card' specifies the number of elements of the 
 ** input set.
 ** The return value is a 32-bit unsigned integer that represents the 
 ** cardinality of the generated power set.
 **/
unsigned powerset(char ***pset, char **set, unsigned card, unsigned k)
{
	char **inset = NULL;
	unsigned pscount = 1;

	while(k-- && pscount)
	{
		if(*pset != NULL)
		{
			if(inset != NULL)
			{
				set_destroy(inset, card);
			}

			inset = *pset;
			card = pscount;
			pscount = powerset_core(pset, inset, card);
		}
		else
		{
			pscount = powerset_core(pset, set, card);
		}
	}

	return pscount;
}

/**
 ** For internal use only. This function contains the logic for generating the
 ** power sets.
 **/
__forceinline unsigned powerset_core(char ***pset, char **set, unsigned card)
{
	char c, *cptr;	
	unsigned pscount = 1 << card, i;
	register unsigned mask;
	size_t buffsize = 0, pi;
	
	*pset = (char**) pool_alloc(pscount * sizeof(char*));

	if(*pset)
	{
		CREATEEMPTYSET((*pset)[0]);
	}
	else
	{
		pscount = 0;
	}

	for(mask = 1; mask < pscount; mask++)
	{
		buffsize = pool_rsv();

		pi = 1;
		c = OPENBRACE;
		pool_write(&c);

		for(i = 0; i < card && pi+2 < buffsize; i++)
		{
			if(mask & (1 << i)) // treat generated integer value as characteristic vector
			{
				if(pi > 1) 
				{
					c = ',';
					pool_write(&c);
					++pi;
				}

				cptr = set[i]; // make a local copy of the pointer!!!

				while(*cptr) // copy element :-D
				{
					pool_write(cptr++);
					++pi;
				}
			}
		}

		if(pi+2 > buffsize)
		{
			pscount = 0;
		}
		else
		{
			c = CLOSEBRACE;
			pool_write(&c);
			c = '\0';
			pool_write(&c);

			(*pset)[mask] = (char*) pool_rls();
		}
	}
	
	return pscount;
}

/**
 ** Writes the specified set to the standard output stream (stdout).
 ** The parameter 'card' specifies the number of elements of the 
 ** input set.
 **/
void set_print(const char **set, const unsigned card)
{
	unsigned i;

	if(card > 0)
	{
		printf("{");

		for(i = 0; i < card; ++i)
		{
			printf(i == (card-1) ? "%s" : "%s,", set[i]);
		}

		printf("}");		
	}
	else
	{
		printf("OUT OF MEMORY :'(");
	}
}

/**
 ** Releases the memory occupied by the specified set.
 ** The freed memory is pushed back to the memory pool.
 **/
void set_destroy(char **set, const unsigned card)
{
	int i = (int)card;

	while(--i >= 0)
	{
		pool_free(set[i]);
		set[i] = NULL;
	}
		
	pool_free(set);
}

/****************** Memory Pool ******************/

typedef union block {
	struct {
		union block* ptr; // pointer to the next free block
		unsigned size; // size of the block (in block units)
	} ci; // contains the control information
	long align; // force 4-byte alignment of blocks
} Block;

static Block pool;
static Block* freep = NULL; // pointer to the first free block
static Block* writep = NULL; // pointer to the block reserved for writing data
static char* writepi; // position indicator used for writing to writep

/**
 ** Creates a memory pool of the specified size bytes of memory.
 ** If the function failed to create the memory pool of the specified size,
 ** the function does nothing.
 **/
void pool_create(const size_t size)
{
	size_t nunits;

	if (pool.ci.ptr) // pool already exists
	{
		return;
	}

	nunits = (size + sizeof(Block) - 1) / sizeof(Block); // round up to the block-sized units
	nunits += 1; // reserve one unit for the control information of the header block

	pool.ci.size = nunits * sizeof(Block);
	pool.ci.ptr = freep = (Block*) malloc(pool.ci.size);
	
	if(freep) // success
	{
		freep->ci.ptr = freep;
		freep->ci.size = nunits;
	}
}

/**
 ** Creates a memory pool utilizing all the available system memory.
 ** This function should not be used for competitive programming.
 **/
void pool_create_maxalloc(void)
{
	register size_t i = 1 << ((sizeof(size_t) << 3) - 1), imin, imax, imid;
	void* ptr = NULL;

	// Use binary search to find the maximum size
	for(imin=(i>>1), imax=i; imin < imax-1;) // ensure the even numbers
	{
		imid = (imin + imax) >> 1; 

		if(ptr = malloc(imid)) 
		{ 
			imin = imid; // search upper subset
			free(ptr); 
		}
		else 
		{
			imax = imid; // search lower subset
		}
	}

	pool_create(imin - 2 * sizeof(Block));
}

/**
 ** Disposes the memory pool and frees the memory previously allocated 
 ** by a call to pool_create.
 **/
void pool_destroy(void)
{
	if(!pool.ci.ptr)
	{
		return;
	}

	free(pool.ci.ptr);
	pool.ci.ptr = freep = writep = NULL;
}

/**
 ** Allocates a block of size bytes of memory, returning a pointer to the 
 ** beginning of the block.
 ** If the function failed to allocate the requested block of memory from the
 ** memory pool, a null pointer is returned.
 **/
void* pool_alloc(const size_t size)
{
	Block *p, *prevp;
	size_t nunits;

	if(!pool.ci.ptr || !freep || writep) return NULL;

	nunits = (size + sizeof(Block) - 1) / sizeof(Block); // round up to the next block size
	nunits += 1; // reserve one unit for the block struct containing the control information

	for(prevp = freep, p = prevp->ci.ptr; ; prevp = p, p = p->ci.ptr)
	{
		if (p->ci.size >= nunits)
		{
			if(p->ci.size == nunits)
			{
				prevp->ci.ptr = p->ci.ptr;
			}
			else // allocate the tail end
			{
				p->ci.size -= nunits;
				p += p->ci.size; // move to the allocated block
				p->ci.size = nunits;
			}

			freep = prevp != p ? prevp : NULL; // adjust the header to the last known free block

			return (void*) (p+1);
		}

		if(p == freep) // out of memory
		{
			return NULL;
		}
	}
}

/**
 ** Frees a block of memory, pointed to by the specified 'mem' pointer,
 ** making it available again for further allocations.
 ** If 'mem' is a null pointer, the function does nothing.
 ** If 'mem' does not point to a block of memory allocated with the 
 ** pool_alloc function, the behavior is undefined.
 **/
void pool_free(void* mem)
{
	Block *bp, *p;

	if(!mem || !pool.ci.ptr || writep) return;

	bp = (Block*) mem - 1; // get a pointer to the block structure

	if(bp->ci.size < 1 || !freep) 
	{
		if(!freep)
		{
			freep = bp;
			freep->ci.ptr = freep;
		}

		return;
	}	

	// Search for the two surrounding blocks
	for(p = freep; bp <= p || bp >= p->ci.ptr; p = p->ci.ptr) 
	{
		if(p >= p->ci.ptr && (bp >= p || bp <= p->ci.ptr))
		{
			break; // deallocate block at the start or end of the list
		}
	}

	// In order to minimize the fragmentation, the adjacent blocks of free memory are
	// combined into a single bigger block.
	if (bp + bp->ci.size == p->ci.ptr) // merge with upper nbr
	{
		bp->ci.size += p->ci.ptr->ci.size;
		bp->ci.ptr = p->ci.ptr->ci.ptr;
	}
	else // link the freed block with the next free block
	{
		bp->ci.ptr = p->ci.ptr; 
	}

	if(p + p->ci.size == bp) // merge with lower nbr
	{
		p->ci.size += bp->ci.size;
		p->ci.ptr = bp->ci.ptr;
	}
	else // link the previous free block with the freed block 
	{
		p->ci.ptr = bp;
	}

	freep = p;
}

/**
 ** Reserves a memory pool for the exclusive write access.
 ** The return value is the size (in bytes) of the largest available 
 ** block of memory that can be used for writing the data.
 **/
size_t pool_rsv(void)
{
	Block *lp, *p;
	size_t avblsize;

	if(!pool.ci.ptr || !freep || writep) return 0;

	for(lp = freep, p = lp->ci.ptr; p != freep; p = p->ci.ptr)
	{
		if(p->ci.size > lp->ci.size)
		{
			lp = p;
		}
	}

	avblsize = (lp->ci.size - 1) * sizeof(Block);

	if(avblsize)
	{
		writep = lp;
		writepi = (char*) (lp + 1);
	}

	return avblsize;
}

/**
 ** Releases a memory pool from the exclusive write access.
 ** The return value is a pointer to the beginning of the block 
 ** allocated by the subsequent calls to the pool_write function.
 **/
void* pool_rls(void)
{
	Block *prevp, *p;
	size_t nunits, funits;

	if(!writep) return NULL;

	// Calculate the size of the block that contains the written data
	nunits = ((writepi - (char*)(writep+1)) + sizeof(Block) - 1) / sizeof(Block) + 1;

	if(nunits-1)
	{
		for(prevp = freep; prevp->ci.ptr != writep; prevp = prevp->ci.ptr);

		if(writep->ci.size == nunits)
		{
			if(prevp == writep)
			{
				freep = NULL;
			}
			else 
			{
				if(writep == freep)
				{
					freep = writep->ci.ptr;
				}

				prevp->ci.ptr = writep->ci.ptr;	
			}		
		}
		else
		{
			p = writep;
			funits = writep->ci.size - nunits;

			writep->ci.size = nunits; // set the size of the newly allocated block

			p += nunits;
			p->ci.size = funits;
			p->ci.ptr = writep->ci.ptr;
			prevp->ci.ptr = p;

			if(writep == freep)
			{
				freep = p;

				if(prevp == writep)
				{
					p->ci.ptr = p;
				}
			}
		}

		p = writep + 1;
	}
	else
	{
		p = NULL;
	}

	writep = NULL;

	return (void*) p;
}

/**
 ** Writes the specified byte into the block of memory reserved by 
 ** the previous call to the pool_rsv function.
 **/
void pool_write(char* c)
{
	if(c && writep)
	{
		*writepi++ = *c;
	}
}
