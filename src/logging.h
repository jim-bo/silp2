/*
 * logging.h
 *
 *  Created on: Feb 28, 2012
 *      Author: eljimbo
 */

#ifndef LOGGING_H_
#define LOGGING_H_

#include <stdio.h>
#include <utility>

// logging and debugging.
#define DEBUG 1
#define VERBOSE 1

static void WRITE_OUT(const char * txt){
	if( VERBOSE == 1 ) {
		fprintf(stdout,"%s",txt);
		fflush(stdout);
	}
}

static void WRITE_ERR(const char * txt){
	fprintf(stderr,"%s",txt);
	fflush(stdout);
}

static bool ACTIVITY(int i, int j, int k){
	if(i != 0 && (i % k) == 0){
		if( VERBOSE == 1 ) fprintf(stdout,"%d of %d\n", i, j);
		fflush(stdout);
		if( DEBUG == 1 ) {
			return true;
		}
	}
	return false;
}


#endif /* LOGGING_H_ */
