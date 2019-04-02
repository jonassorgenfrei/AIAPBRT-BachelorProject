#include <stdio.h>
#include <vector>

/// <summary>
/// Mains the specified argc.
/// </summary>
/// <param name="argc">The argc.</param>
/// <param name="argv">The argv.</param>
/// <returns></returns>
int main(int argc, char *argv[]) {
	Options options;
	std::vector<std::string> filenames;
	
	// Process command-line arguments
	for (int i = 1; i < argc; ++i) {
		if (!strcmp(argv[i], "--ncores") || !strcmp(argv[i], "--nthreads"))
			options.nThreads = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--outfile"))
			options.quickRender = true;
		else if (!strcmp(argv[i], "--quiet"))
			options.quiet = true;
		else if (!strcmp(argv[i], "--verbose"))
			options.verbose = true;
		else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")){
			printf("usage: pbrt [--nthreads n] [--outfile filename] [--quick] [--quiet] "
				   "[--verbose] [--help] <filename.pbrt> ...\n");
			return 0;
		}
		else filenames.push_back(argv[i]);
	}

	pbrtInit(options);
	
	if (filenames.size() == 0) {
		// Parse Scene from standard input
		ParseFile("-");
	} else {
		/* Parse scene from input files */
		for(const std::string &f : filenames)
			if (!ParseFile(f)) 
				Error("Couldn't open scene file \"%s\"", f.c_str());
			;
	}

	//pbrtCleanup();
	return 0;
}