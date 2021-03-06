/*
This file is intended to serve as a usage example.
Call with:
./main [filename]
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>


#include "Graph.h"
#include "graph_simple.h"
#include "read_table.h"

//for time measurement
#define dt(b,a) ( (double) (a.tv_usec - b.tv_usec) / 1000000 + (double) (a.tv_sec - b.tv_sec))


int main(int argc, char** argv){
	char* fn = 0;
	bool read_part = false;
	for(int i=1;i<argc;i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'i':
			fn = argv[i+1];
			if(fn[0] == '-' && fn[1] == 0) fn = 0;
			i++;
			break;
		case 'p':
			read_part = true;
			break;
		default:
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
	}
	FILE* fp;
	if(fn) fp = fopen(fn,"r");
	else fp = stdin;
	graph g1;
	std::vector<unsigned int> ids;
	if(g1.read_graph(fp,read_part,&ids)) {
		fprintf(stderr,"Error reading the graph!\n");
		return 1;
	}
	if(g1.make_symmetric()) {
		fprintf(stderr,"Error creating a symmetric graph!\n");
		return 1;
	}
	MVGraph g(g1);
	
	//~ fprintf(stderr,"%u\n",g.matchnum);
	
	struct timeval start,end;  
	gettimeofday(&start, NULL);
	g.greedy_init();
	g.max_match();
	gettimeofday(&end, NULL);
	
	fprintf(stderr,"%f\n",dt(start,end));
	fprintf(stderr,"%u\n",g.matchnum);
	
	g.write_matches(stdout,[ids](nodeid i){ return ids[i]; });
	
	return 0;
}
