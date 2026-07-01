
class CUDATimer {
private:
	cudaEvent_t start_, stop_;

public:
	CUDATimer() {
		cudaEventCreate(&start_);
		cudaEventCreate(&stop_);
		cudaEventRecord(start_, 0);
	}

	float stop() {
		cudaEventRecord(stop_, 0);
		cudaEventSynchronize(stop_);
		float time;
		cudaEventElapsedTime(&time, start_, stop_);
		return time * 1e-3;
	}
};

class CPUTimer {
private:
	double start;

	double getTime() {
		timespec clock;
		clock_gettime(CLOCK_MONOTONIC, &clock);
		return clock.tv_sec + clock.tv_nsec * 1e-9;
	}

public:
	CPUTimer() : start(getTime()) { }

	double stop() {
		return getTime() - start;
	}
};
