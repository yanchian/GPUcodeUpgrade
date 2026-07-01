
#define checkForError() \
  do { \
  	cudaError_t error = cudaGetLastError(); \
	  if (error != cudaSuccess) { \
  		std::cerr << "CUDA error " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString(error) << std::endl; \
	  	exit(1); \
    } \
  } while (false);

