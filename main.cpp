#include "Utility.h"
#include "Graph.h"

void print_usage() {
	printf("Usage: [1]exe [2]graph-dir\n");
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		print_usage();
		return 0;
	}


#ifdef _LINUX_
	struct timeval start, end1, end;
	gettimeofday(&start, NULL);
#endif
	
	Graph *graph = new Graph(argv[1]);
	graph->read_graph_GIS();
	// graph->read_graph_GIS2();
#ifdef _LINUX_
	gettimeofday(&end1, NULL);

	long long mtime1, seconds1, useconds1;
	seconds1 = end1.tv_sec - start.tv_sec;
	useconds1 = end1.tv_usec - start.tv_usec;
	mtime1 = seconds1*1000000 + useconds1;
#endif

	graph->GIS();
#ifdef _LINUX_
	gettimeofday(&end, NULL);

	long long mtime, seconds, useconds;
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = seconds*1000000 + useconds;

	//printf("Total time excluding IO is: %lld\n", mtime-mtime1);
#endif

	return 0;
}
