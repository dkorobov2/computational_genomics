#pragma warning(disable:4996)
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <CL/cl.h>
//#include <sys/mman.h>
#include "data_structures.h"
#include "params.h"
#include <cstdlib>
#include <math.h>

using namespace std;

void checkError(cl_int err)
{
	if (err != CL_SUCCESS) {
		printf("Error: %d\n", err);
		system("pause");
	}
}

/* Create program from a file and compile it */
cl_program build_program(cl_context ctx, cl_device_id dev, cl_uint  num_of_devices, const char* filename) {

	cl_program program;
	FILE *program_handle;
	char *program_buffer, *program_log;
	size_t program_size, log_size;
	int err;

	/* Read program file and place content into buffer */
	program_handle = fopen(filename, "rb");
	if (program_handle == NULL) {
		perror("Couldn't find the program file");
		exit(1);
	}
	fseek(program_handle, 0, SEEK_END);
	program_size = ftell(program_handle);
	rewind(program_handle);

	program_buffer = (char*)malloc(program_size + 1);

	fread(program_buffer, sizeof(char), program_size, program_handle);
	//printf("%d\n", program_size);
	//system("pause");
	program_buffer[program_size] = '\0';


	//printf("%s\n", program_buffer);
	fclose(program_handle);

	/* Create program from file */
	program = clCreateProgramWithSource(ctx, 1, (const char **)&program_buffer, NULL, &err);
	if (err < 0) {
		perror("Couldn't create the program");
		exit(1);
	}
	free(program_buffer);

	/* Build program */
	//err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

	err = clBuildProgram(program, num_of_devices, &dev, NULL, NULL, NULL);

	if (err < 0) {

		/* Find size of log and print to std output */
		clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
		program_log = (char*)malloc(log_size + 1);
		program_log[log_size] = '\0';
		clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG, log_size + 1, program_log, NULL);
		printf("%s\n", program_log);
		system("pause");
		free(program_log);
		system("pause");
		exit(1);
	}

	return program;
}

void runEntropyKernel(cl_context context, cl_command_queue command_queue, cl_program program, cl_uint num_of_devices, cl_device_id device_id, cl_kernel kernel, size_t global, size_t local, int case_num, ncount_t * counts, float * entropy_out, cl_int grain_size, cl_mem &counts_dev, cl_mem &entropy_out_dev)
{
	cl_int err;
	//printf("%d", sizeof(reference));
	if (case_num == 1)
		clEnqueueWriteBuffer(command_queue, counts_dev, CL_TRUE, 0, (NUM_NUCL + LOCAL) * sizeof(ncount_t), counts, 0, NULL, NULL);
	else if (case_num == 2)
		clEnqueueWriteBuffer(command_queue, counts_dev, CL_TRUE, 0, (NUM_NUCL + 2 * LOCAL) * sizeof(ncount_t), counts, 0, NULL, NULL);
	else if (case_num == 3)
		clEnqueueWriteBuffer(command_queue, counts_dev, CL_TRUE, 0, (NUM_NUCL + LOCAL) * sizeof(ncount_t), counts, 0, NULL, NULL);
	
	program = build_program(context, device_id, num_of_devices, ENTROPY);

	kernel = clCreateKernel(program, "calculateEntropy", &err);
	checkError(err);

	// set the argument list for the kernel command
	err |= clSetKernelArg(kernel, 0, sizeof(cl_mem), &entropy_out_dev);
	err |= clSetKernelArg(kernel, 1, sizeof(unsigned int), &counts_dev);
	err |= clSetKernelArg(kernel, 2, sizeof(cl_int), &grain_size);
	checkError(err);

	cout << "global: :" << global << "local: " << local << endl;
	err = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global, &local, 0, NULL, NULL);
	checkError(err);

	clFinish(command_queue);

	clEnqueueReadBuffer(command_queue, entropy_out_dev, CL_TRUE, 0, NUM_NUCL * sizeof(float), entropy_out, 0, NULL, NULL);

	printEntropy(entropy_out, 10);
}

void runConvolutionKernel(cl_context context, cl_command_queue command_queue, cl_program program2, cl_uint num_of_devices, cl_device_id device_id, cl_kernel kernel, size_t global, size_t local, float * filter, float * convolution_out, cl_uint filter_length, cl_uint input_length, cl_mem &filter_dev, cl_mem &convolution_out_dev, cl_mem &entropy_out_dev)
{
	cl_int err;

	program2 = build_program(context, device_id, num_of_devices, CONVOLUTION);

	kernel = clCreateKernel(program2, "Convolve", &err);
	checkError(err);

	// set the argument list for the kernel command
	err |= clSetKernelArg(kernel, 0, sizeof(cl_mem), &entropy_out_dev);
	err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &filter_dev);
	err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &convolution_out_dev);
	err |= clSetKernelArg(kernel, 3, sizeof(cl_uint), &input_length);
	err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &filter_length);
	checkError(err);

	err = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global, &local, 0, NULL, NULL);
	checkError(err);

	clFinish(command_queue);

	clEnqueueReadBuffer(command_queue, convolution_out_dev, CL_TRUE, 0, NUM_NUCL * sizeof(float), convolution_out, 0, NULL, NULL);

	printConvolution(convolution_out, 10);
}

void runThresholdKernel(cl_context context, cl_command_queue command_queue, cl_program program3, cl_uint num_of_devices, cl_device_id device_id, cl_kernel kernel, size_t global, size_t local, int iteration, unsigned char * threshold_out, float threshold, cl_mem &threshold_out_dev, cl_mem &convolution_out_dev)
{
	cl_int err;

	program3 = build_program(context, device_id, num_of_devices, THRESHOLD);

	kernel = clCreateKernel(program3, "threshold", &err);
	checkError(err);

	// set the argument list for the kernel command
	err |= clSetKernelArg(kernel, 0, sizeof(cl_mem), &convolution_out_dev);
	err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &threshold_out_dev);
	err |= clSetKernelArg(kernel, 2, sizeof(float), &threshold);

	err = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global, &local, 0, NULL, NULL);
	checkError(err);

	clFinish(command_queue);
	
	//first iteration
	if (iteration == 0)
		clEnqueueReadBuffer(command_queue, threshold_out_dev, CL_TRUE, 0, NUM_NUCL * sizeof(unsigned char), threshold_out, 0, NULL, NULL);
	// some iteration in the middle
	// ignore first 512 elements since they were used just used for convolution
	else if (iteration > 0 && iteration < NUM_ITERATIONS - 1)
		clEnqueueReadBuffer(command_queue, threshold_out_dev, CL_TRUE, LOCAL, NUM_NUCL * sizeof(unsigned char), threshold_out, 0, NULL, NULL);
	// final iteration
	// ignore first 512 elements since they were used just used for convolution
	else
		clEnqueueReadBuffer(command_queue, threshold_out_dev, CL_TRUE, LOCAL, NUM_NUCL * sizeof(unsigned char), threshold_out, 0, NULL, NULL);
	
	printThreshold(threshold_out, 10);
}

int main()
{
	printf("Testing Entropy/Convolution\n");
	cl_context context;
	cl_context_properties properties[3];
	cl_kernel kernel = NULL;
	cl_command_queue command_queue;
	cl_program program = NULL, program2 = NULL, program3 = NULL;
	cl_int err;
	cl_uint num_of_platforms = 0;
	cl_platform_id platform_id;
	cl_device_id device_id;
	cl_uint num_of_devices = 0;
	cl_mem counts_dev = NULL, entropy_out_dev = NULL, filter_dev = NULL, convolution_out_dev = NULL, threshold_out_dev = NULL;

	// initialize the Rtree used for storing the start and end positions of each short read
	SRTree sr_tree;

	// initialize arrays from alignment
	uint32_t * nucs = createNucs(NUM_READS);
	uint32_t * cigars = createCigars(NUM_READS);
	uint32_t * starting = createStarting(NUM_READS);

	vector<ncount_t> freq(NUM_NUCL_TOTAL);
	vector<unsigned char> threshold_vector(NUM_NUCL_TOTAL);
	float * entropy_out = new float[NUM_NUCL_TOTAL];
	float * convolution_out = new float[NUM_NUCL_TOTAL];
	//unsigned char * threshold_out = new unsigned char[NUM_NUCL_TOTAL];
	float filter[] = { .1, .2, .4, .2, .1 };

	// mmap file
/*
	fseek(fp, 0L, SEEK_END);
	int file_size = ftell(fp);
	mmap(NULL, file_size, PROT_READ, MAP_ANONYMOUS, fd, 0);
	*/

	// Generate RTree, sr_indexes and counts array
	vector<uint32_t> sr_indexes(NUM_NUCL_TOTAL);
	createRTreeAndFrequencies(nucs, starting, cigars, sr_tree, freq, sr_indexes);

	// Have a pointer to the vectors for kernels
	ncount_t * counts = &freq[0];
	unsigned char * threshold_out = &threshold_vector[0];

	printCounts(counts, 10);
	size_t global, local;

	// retreive a list of platforms avaible
	if (clGetPlatformIDs(1, &platform_id, &num_of_platforms) != CL_SUCCESS)
	{
		printf("Unable to get platform_id\n");
		system("pause");
		return 1;
	}
	printf("platform_id: %i \n", platform_id);


	if (clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &num_of_devices) != CL_SUCCESS)
	{
		printf("Unable to get device_id\n");
		system("pause");
		return 1;
	}

	// context properties list - must be terminated with 0
	properties[0] = CL_CONTEXT_PLATFORM;
	properties[1] = (cl_context_properties)platform_id;
	properties[2] = 0;

	// create a context with the GPU device
	context = clCreateContext(properties, 1, &device_id, NULL, NULL, &err);

	// create command queue using the context and device
	command_queue = clCreateCommandQueue(context, device_id, 0, &err);

	///////////////// BUFFER INITIALIZATION ////////////////////
	// create all the necessary opencl buffers
	counts_dev = clCreateBuffer(context, CL_MEM_READ_ONLY, (NUM_NUCL + 2 * LOCAL)* sizeof(ncount_t), NULL, &err);
	checkError(err);
	entropy_out_dev = clCreateBuffer(context, CL_MEM_WRITE_ONLY, (NUM_NUCL + 2 * LOCAL)* sizeof(float), NULL, &err);
	checkError(err);

	filter_dev = clCreateBuffer(context, CL_MEM_READ_ONLY, FILTER_LENGTH * 4, NULL, &err);
	checkError(err);
	clEnqueueWriteBuffer(command_queue, filter_dev, CL_TRUE, 0, FILTER_LENGTH * sizeof(float), filter, 0, NULL, NULL);
	
	convolution_out_dev = clCreateBuffer(context, CL_MEM_WRITE_ONLY, (NUM_NUCL + 2 * LOCAL) * sizeof(float), NULL, &err);
	checkError(err);

	threshold_out_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, (NUM_NUCL + 2 * LOCAL), NULL, &err);
	checkError(err);

	//begin iterations
	for (int i = 0; i < NUM_ITERATIONS; i++) {
		/////////////////////////////// ENTROPY /////////////////////////////////
		cl_int grain_size = GRAIN_SIZE;
		local = LOCAL;

		// if this is the first batch of elements
		if (i == 0) {
			global = (NUM_NUCL + LOCAL) / GRAIN_SIZE;
			runEntropyKernel(context, command_queue, program, num_of_devices, device_id, kernel, global, local, 1, &counts[i * NUM_NUCL], &entropy_out[i * NUM_NUCL], grain_size, counts_dev, entropy_out_dev);
		}
		// if this is from the middle batch of elements
		else if (i > 0 && i < NUM_ITERATIONS - 1) {
			global = (NUM_NUCL + LOCAL * 2) / GRAIN_SIZE;
			runEntropyKernel(context, command_queue, program, num_of_devices, device_id, kernel, global, local, 2, &counts[i * NUM_NUCL - LOCAL], &entropy_out[i * NUM_NUCL - LOCAL], grain_size, counts_dev, entropy_out_dev);
		}
		// if this is the last batch of elements
		else {
			global = (NUM_NUCL + LOCAL) / GRAIN_SIZE;
			runEntropyKernel(context, command_queue, program, num_of_devices, device_id, kernel, global, local, 3, &counts[i * NUM_NUCL - LOCAL], &entropy_out[i * NUM_NUCL - LOCAL], grain_size, counts_dev, entropy_out_dev);
		}
		/////////////////////////// CONVOLUTION /////////////////////////////
		cl_uint filter_length = FILTER_LENGTH;
		cl_uint input_length = NUM_NUCL;
		local = LOCAL;

		if (i == 0) {
			global = (NUM_NUCL + LOCAL);
			runConvolutionKernel(context, command_queue, program2, num_of_devices, device_id, kernel, global, local, filter, &convolution_out[i * NUM_NUCL], filter_length, input_length, filter_dev, convolution_out_dev, entropy_out_dev);
		}
		else if (i > 0 && i < NUM_ITERATIONS - 1) {
			global = (NUM_NUCL + LOCAL * 2);
			runConvolutionKernel(context, command_queue, program2, num_of_devices, device_id, kernel, global, local, filter, &convolution_out[i * NUM_NUCL - LOCAL], filter_length, input_length, filter_dev, convolution_out_dev, entropy_out_dev);
		}
		else {
			global = (NUM_NUCL + LOCAL);
			runConvolutionKernel(context, command_queue, program2, num_of_devices, device_id, kernel, global, local, filter, &convolution_out[i * NUM_NUCL - LOCAL], filter_length, input_length, filter_dev, convolution_out_dev, entropy_out_dev);
		}
		/////////////////// CPU CONVOLUTION CHECK //////////////////
	/*
		float * convolution_out_cpu = new float[NUM_NUCL];
		convolve(entropy_out, NUM_NUCL, filter, FILTER_LENGTH, convolution_out_cpu);
		printConvolution(convolution_out_cpu, 10);
	*/

	//////////////////////// THRESHOLD /////////////////////////
		float threshold = THRESHOLD_VAL;
		global = NUM_NUCL;
		local = LOCAL;

		if (i == 0) {
			global = (NUM_NUCL + LOCAL);
			runThresholdKernel(context, command_queue, program3, num_of_devices, device_id, kernel, global, local, i, &threshold_out[i * NUM_NUCL], threshold, threshold_out_dev, convolution_out_dev);
		}
		else if (i > 0 && i < NUM_ITERATIONS - 1) {
			global = (NUM_NUCL + LOCAL * 2);
			runThresholdKernel(context, command_queue, program3, num_of_devices, device_id, kernel, global, local, i, &threshold_out[i * NUM_NUCL], threshold, threshold_out_dev, convolution_out_dev);
		}
		else {
			global = (NUM_NUCL + LOCAL);
			runThresholdKernel(context, command_queue, program3, num_of_devices, device_id, kernel, global, local, i, &threshold_out[i * NUM_NUCL], threshold, threshold_out_dev, convolution_out_dev);
		}
	}
	
	// done with kernels, create the de bruijn input
	vector<uint32_t> DBInput = createDBInput(nucs, threshold_vector, sr_indexes, sr_tree);
	
	// print a sample of the debruin input
	print_elements(DBInput, 9);

	delete[] nucs;
	delete[] cigars;
	delete[] starting;
	delete []entropy_out;
	delete []convolution_out;

	clReleaseMemObject(counts_dev);
	clReleaseMemObject(filter_dev);
	clReleaseMemObject(convolution_out_dev);
	clReleaseMemObject(threshold_out_dev);

	clReleaseProgram(program);
	clReleaseProgram(program2);
	clReleaseProgram(program3);
	clReleaseKernel(kernel);
	system("pause");
	return 0;
}