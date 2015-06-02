static char adSid[]="$Id: adStack.c,v 1.4 10/.0/.0 .1:.5:.2 llh Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define ONE_BLOCK_SIZE 16384
#define CHUNK_SIZE 4096

/* The main stack is a double-chain of DoubleChainedBlock objects.
 * Each DoubleChainedBlock holds an array[ONE_BLOCK_SIZE] of char. */
typedef struct _doubleChainedBlock{
  struct _doubleChainedBlock *prev ;
  char                       *stackBottom ;
  char                       *stackTop ;
  char                       *stackLimit ;
  struct _doubleChainedBlock *next ;
} DoubleChainedBlock ;

/* Globals that define the current position in the stack: */
static DoubleChainedBlock *curStack = NULL ;
static char               *curStackBottom = NULL ;
static char               *curStackTop    = NULL ;
static char               *curStackLimit  = NULL ;
/* Globals that define the current LOOKing position in the stack: */
static DoubleChainedBlock *lookStack = NULL ;
static char               *lookStackBottom = NULL ;
static char               *lookStackTop    = NULL ;
static char               *lookStackLimit  = NULL ;

/* Used before PUSHing.
 * Resets the LOOKing position if it was active.
 * Checks that there is enough space left to hold "nbChars" chars.
 * Otherwise, allocates the necessary space. */
void check(int nbChars) {
  if (lookStack) lookStack = NULL ;
  if (curStackTop+nbChars > curStackLimit) {
    if ((curStack == NULL) || (curStack->next == NULL)) {
      DoubleChainedBlock *newStack ;
      char *contents = (char*)malloc(ONE_BLOCK_SIZE*sizeof(char)) ;
      newStack = (DoubleChainedBlock*)malloc(sizeof(DoubleChainedBlock)) ;
      if ((contents == NULL) || (newStack == NULL)) {
	  DoubleChainedBlock *stack = curStack ;
	  int nbBlocks = 0 ;
	  while(stack) {
	      stack = stack->prev ;
	      nbBlocks++ ;
	  }
	  printf("Out of memory (allocated %i blocks of %i bytes)\n",
		 nbBlocks, ONE_BLOCK_SIZE) ;
          exit(0);
      }
      newStack->prev = curStack ;
      if (curStack != NULL) {
	curStack->stackTop = curStackTop ;
	curStack->next = newStack ;
      }
      newStack->next = NULL ;
      newStack->stackBottom = contents ;
      newStack->stackTop = contents ;
      newStack->stackLimit = contents+ONE_BLOCK_SIZE ;
      curStack = newStack ;
    } else {
      curStack->stackTop = curStackTop ;
      curStack = curStack->next ;
    }
    curStackBottom = curStack->stackBottom ;
    curStackTop    = curStackBottom ;
    curStackLimit  = curStack->stackLimit ;
  }
}
/* Used before POPping.
 * Resets the LOOKing position if it was active.
 * Checks that the current block is not empty.
 * If it is, put the pointer back to the previous block. */
void checkBack() {
  if (lookStack) lookStack = NULL ;
  if (curStackTop == curStackBottom) {
    curStack->stackTop = curStackBottom ;
    curStack = curStack->prev ;
    curStackBottom = curStack->stackBottom ;
    curStackTop    = curStack->stackTop ;
    curStackLimit  = curStack->stackLimit ;
  }
}
/* Used before LOOKing.
 * Activates the LOOKing position if it was reset.
 * Checks that the current LOOKing block is not empty.
 * If it is, put the LOOK pointer back to the previous block. */
void checkLookBack() {
  if (lookStack == NULL) {
    lookStack = curStack ;
    lookStackBottom = curStackBottom ;
    lookStackTop = curStackTop ;
    lookStackLimit = curStackLimit ;
  }
  if (lookStackTop == lookStackBottom) {
    lookStack->stackTop = lookStackBottom ;
    lookStack = lookStack->prev ;
    lookStackBottom = lookStack->stackBottom ;
    lookStackTop    = lookStack->stackTop ;
    lookStackLimit  = lookStack->stackLimit ;
  }
}

/* PUSHes "nbChars" consecutive chars,
 * from a location starting at address "x".
 * nbChars is assumed no larger than CHUNK_SIZE */
void pushN(char *x, int nbChars) {
  check(nbChars) ;
  memcpy(curStackTop,x,nbChars);
  curStackTop+=nbChars ;
}
/* POPs "nbChars" consecutive chars,
 * to a location starting at address "x".
 * nbChars is assumed no larger than CHUNK_SIZE */
void popN(char *x, int nbChars) {
  checkBack() ;
  curStackTop-=nbChars ;
  memcpy(x,curStackTop,nbChars);
}
/* LOOKs "nbChars" consecutive chars,
 * to a location starting at address "x".
 * LOOKing is just like POPping, except that the main pointer
 * remains in place, so that the value is not POPped.
 * Further PUSHs or POPs will start from the same place as if
 * no LOOK had been made.
 * nbChars is assumed no larger than CHUNK_SIZE */
void lookN(char *x, int nbChars) {
  checkLookBack() ;
  lookStackTop-=nbChars ;
  memcpy(x,lookStackTop,nbChars);
}

/* PUSHes a large number "n" of consecutive chars,
 * from a location starting at address "x".
 * This "n"-sized array is cut into pieces no larger than CHUNK_SIZE */
void pushNarray(char *x, int n) {
  int tailSize = n%CHUNK_SIZE ;
  char *xmax = x+n-tailSize ;
  char *xin = x ;
  while(xin<xmax) {
    pushN(xin, CHUNK_SIZE) ;
    xin += CHUNK_SIZE ;
  }
  if (tailSize>0) pushN(xin, tailSize) ;
}
/* POPs a large number "n" of consecutive chars,
 * to a location starting at address "x".
 * This "n"-sized array is cut into pieces no larger than CHUNK_SIZE */
void popNarray(char *x, int n) {
  int tailSize = n%CHUNK_SIZE ;
  char *xin = x+n-tailSize ;
  if (tailSize>0) popN(xin, tailSize) ;
  while(xin>x) {
    xin -= CHUNK_SIZE ;
    popN(xin, CHUNK_SIZE) ;
  }
}
/* LOOKs a large number "n" of consecutive chars,
 * to a location starting at address "x".
 * This "n"-sized array is cut into pieces no larger than CHUNK_SIZE */
void lookNarray(char *x, int n) {
  int tailSize = n%CHUNK_SIZE ;
  char *xin = x+n-tailSize ;
  if (tailSize>0) lookN(xin, tailSize) ;
  while(xin>x) {
    xin -= CHUNK_SIZE ;
    lookN(xin, CHUNK_SIZE) ;
  }
}

/********** Exported PUSH/POP/LOOK functions : ************/

void pushcharacter_(char *x) {
  pushN(x,1) ;
}
void popcharacter_(char *x) {
  popN(x,1) ;
}
void lookcharacter_(char *x) {
  lookN(x,1) ;
}

void pushboolean_(char *x) {
  pushN(x,4) ;
}
void popboolean_(char *x) {
  popN(x,4) ;
}
void lookboolean_(char *x) {
  lookN(x,4) ;
}

void pushinteger4_(char *x) {
  pushN(x,4) ;
}
void popinteger4_(char *x) {
  popN(x,4) ;
}
void lookinteger4_(char *x) {
  lookN(x,4) ;
}

void pushinteger8_(char *x) {
  pushN(x,8) ;
}
void popinteger8_(char *x) {
  popN(x,8) ;
}
void lookinteger8_(char *x) {
  lookN(x,8) ;
}

void pushinteger16_(char *x) {
  pushN(x,16) ;
}
void popinteger16_(char *x) {
  popN(x,16) ;
}
void lookinteger16_(char *x) {
  lookN(x,16) ;
}

void pushreal4_(char *x) {
  pushN(x,4) ;
}
void popreal4_(char *x) {
  popN(x,4) ;
}
void lookreal4_(char *x) {
  lookN(x,4) ;
}

void pushreal8_(char *x) {
  pushN(x,8) ;
}
void popreal8_(char *x) {
  popN(x,8) ;
}
void lookreal8_(char *x) {
  lookN(x,8) ;
}

void pushreal16_(char *x) {
  pushN(x,16) ;
}
void popreal16_(char *x) {
  popN(x,16) ;
}
void lookreal16_(char *x) {
  lookN(x,16) ;
}

void pushreal32_(char *x) {
  pushN(x,32) ;
}
void popreal32_(char *x) {
  popN(x,32) ;
}
void lookreal32_(char *x) {
  lookN(x,32) ;
}

void pushcomplex4_(char *x) {
  pushN(x,4) ;
}
void popcomplex4_(char *x) {
  popN(x,4) ;
}
void lookcomplex4_(char *x) {
  lookN(x,4) ;
}

void pushcomplex8_(char *x) {
  pushN(x,8) ;
}
void popcomplex8_(char *x) {
  popN(x,8) ;
}
void lookcomplex8_(char *x) {
  lookN(x,8) ;
}

void pushcomplex16_(char *x) {
  pushN(x,16) ;
}
void popcomplex16_(char *x) {
  popN(x,16) ;
}
void lookcomplex16_(char *x) {
  lookN(x,16) ;
}

void pushcomplex32_(char *x) {
  pushN(x,32) ;
}
void popcomplex32_(char *x) {
  popN(x,32) ;
}
void lookcomplex32_(char *x) {
  lookN(x,32) ;
}

/******************** The same for arrays: ****************/

void pushcharacterarray_(char *x, int *n) {
  pushNarray(x,*n) ;
}
void popcharacterarray_(char *x, int *n) {
  popNarray(x,*n) ;
}
void lookcharacterarray_(char *x, int *n) {
  lookNarray(x,*n) ;
}

void pushbooleanarray_(char *x, int *n) {
  pushNarray(x,(*n*4)) ;
}
void popbooleanarray_(char *x, int *n) {
  popNarray(x,(*n*4)) ;
}
void lookbooleanarray_(char *x, int *n) {
  lookNarray(x,(*n*4)) ;
}

void pushinteger4array_(char *x, int *n) {
  pushNarray(x,(*n*4)) ;
}
void popinteger4array_(char *x, int *n) {
  popNarray(x,(*n*4)) ;
}
void lookinteger4array_(char *x, int *n) {
  lookNarray(x,(*n*4)) ;
}

void pushinteger8array_(char *x, int *n) {
  pushNarray(x,(*n*8)) ;
}
void popinteger8array_(char *x, int *n) {
  popNarray(x,(*n*8)) ;
}
void lookinteger8array_(char *x, int *n) {
  lookNarray(x,(*n*8)) ;
}

void pushinteger16array_(char *x, int *n) {
  pushNarray(x,(*n*16)) ;
}
void popinteger16array_(char *x, int *n) {
  popNarray(x,(*n*16)) ;
}
void lookinteger16array_(char *x, int *n) {
  lookNarray(x,(*n*16)) ;
}

void pushreal4array_(char *x, int *n) {
  pushNarray(x,(*n*4)) ;
}
void popreal4array_(char *x, int *n) {
  popNarray(x,(*n*4)) ;
}
void lookreal4array_(char *x, int *n) {
  lookNarray(x,(*n*4)) ;
}

void pushreal8array_(char *x, int *n) {
  pushNarray(x,(*n*8)) ;
}
void popreal8array_(char *x, int *n) {
  popNarray(x,(*n*8)) ;
}
void lookreal8array_(char *x, int *n) {
  lookNarray(x,(*n*8)) ;
}

void pushreal16array_(char *x, int *n) {
  pushNarray(x,(*n*16)) ;
}
void popreal16array_(char *x, int *n) {
  popNarray(x,(*n*16)) ;
}
void lookreal16array_(char *x, int *n) {
  lookNarray(x,(*n*16)) ;
}

void pushreal32array_(char *x, int *n) {
  pushNarray(x,(*n*32)) ;
}
void popreal32array_(char *x, int *n) {
  popNarray(x,(*n*32)) ;
}
void lookreal32array_(char *x, int *n) {
  lookNarray(x,(*n*32)) ;
}

void pushcomplex4array_(char *x, int *n) {
  pushNarray(x,(*n*4)) ;
}
void popcomplex4array_(char *x, int *n) {
  popNarray(x,(*n*4)) ;
}
void lookcomplex4array_(char *x, int *n) {
  lookNarray(x,(*n*4)) ;
}

void pushcomplex8array_(char *x, int *n) {
  pushNarray(x,(*n*8)) ;
}
void popcomplex8array_(char *x, int *n) {
  popNarray(x,(*n*8)) ;
}
void lookcomplex8array_(char *x, int *n) {
  lookNarray(x,(*n*8)) ;
}

void pushcomplex16array_(char *x, int *n) {
  pushNarray(x,(*n*16)) ;
}
void popcomplex16array_(char *x, int *n) {
  popNarray(x,(*n*16)) ;
}
void lookcomplex16array_(char *x, int *n) {
  lookNarray(x,(*n*16)) ;
}

void pushcomplex32array_(char *x, int *n) {
  pushNarray(x,(*n*32)) ;
}
void popcomplex32array_(char *x, int *n) {
  popNarray(x,(*n*32)) ;
}
void lookcomplex32array_(char *x, int *n) {
  lookNarray(x,(*n*32)) ;
}

/************* Debug displays of the state of the stack: ***********/

void printtopplace_() {
    DoubleChainedBlock *stack = curStack ;
    int nbBlocks = 0 ;
    while(stack) {
	stack = stack->prev ;
	nbBlocks++ ;
    }
    printf("Stack  top: %i+%i\n",nbBlocks,curStackTop - curStackBottom) ;
}

void printlookingplace_() {
    if (lookStack == NULL)
	printtopplace_() ;
    else {
	DoubleChainedBlock *stack = lookStack ;
	int nbBlocks = 0 ;
	while(stack) {
	    stack = stack->prev ;
	    nbBlocks++ ;
	}
	printf("Stack look: %i+%i\n",nbBlocks,lookStackTop - lookStackBottom) ;
    }
}
