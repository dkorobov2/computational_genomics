#pragma warning(disable:4996)
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <CL/cl.h>
#include "data_structures.h"
#include "params.h"
#include <math.h>
#include <chrono>

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
	printf("%d\n", program_size);
	//system("pause");
	program_buffer[program_size] = '\0';


	printf("%s\n", program_buffer);
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

void runEntropyKernel(cl_context context, cl_command_queue command_queue, cl_program program, cl_uint num_of_devices, cl_device_id device_id, cl_kernel kernel, size_t global, size_t local, ncount_t * counts, float * entropy_out, cl_int grain_size, cl_mem &counts_dev, cl_mem &entropy_out_dev)
{
	cl_int err = 0;
/*	counts_dev = clCreateBuffer(context, CL_MEM_READ_ONLY, NUM_NUCL * sizeof(ncount_t), NULL, &err);
	checkError(err);
	entropy_out_dev = clCreateBuffer(context, CL_MEM_WRITE_ONLY, NUM_NUCL * sizeof(float), NULL, &err);
	checkError(err);
	//printf("%d", sizeof(reference));
	clEnqueueWriteBuffer(command_queue, counts_dev, CL_TRUE, 0, NUM_NUCL * sizeof(ncount_t), counts, 0, NULL, NULL);

	program = build_program(context, device_id, num_of_devices, ENTROPY);

	kernel = clCreateKernel(program, "calculateEntropy", &err);
	checkError(err);
*/
	// set the argument list for the kernel command
	err |= clSetKernelArg(kernel, 0, sizeof(cl_mem), &entropy_out_dev);
	checkError(err);
	err |= clSetKernelArg(kernel, 1, sizeof(unsigned int), &counts_dev);
	checkError(err);
	err |= clSetKernelArg(kernel, 2, sizeof(cl_int), &grain_size);
	checkError(err);

	err = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global, &local, 0, NULL, NULL);
	checkError(err);

	clFinish(command_queue);

	//clEnqueueReadBuffer(command_queue, entropy_out_dev, CL_TRUE, 0, NUM_NUCL * sizeof(float), entropy_out, 0, NULL, NULL);

	//printEntropy(entropy_out, 10);
}

void runConvolutionKernel(cl_context context, cl_command_queue command_queue, cl_program program2, cl_uint num_of_devices, cl_device_id device_id, cl_kernel kernel, size_t global, size_t local, float * filter, float * convolution_out, cl_uint filter_length, cl_uint input_length, cl_mem &filter_dev, cl_mem &convolution_out_dev, cl_mem &entropy_out_dev)
{
	cl_int err = 0;
	filter_dev = clCreateBuffer(context, CL_MEM_READ_ONLY, FILTER_LENGTH * 4, NULL, &err);
	checkError(err);
	convolution_out_dev = clCreateBuffer(context, CL_MEM_WRITE_ONLY, NUM_NUCL * sizeof(float), NULL, &err);
	checkError(err);

	clEnqueueWriteBuffer(command_queue, filter_dev, CL_TRUE, 0, FILTER_LENGTH * sizeof(float), filter, 0, NULL, NULL);

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

void runThresholdKernel(cl_context context, cl_command_queue command_queue, cl_program program3, cl_uint num_of_devices, cl_device_id device_id, cl_kernel kernel, size_t global, size_t local, cl_uint * threshold_out, float threshold, cl_mem &threshold_out_dev, cl_mem &convolution_out_dev)
{
	cl_int err = 0;
	threshold_out_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, NUM_NUCL * 4, NULL, &err);
	checkError(err);

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

	clEnqueueReadBuffer(command_queue, threshold_out_dev, CL_TRUE, 0, NUM_NUCL * sizeof(cl_uint), threshold_out, 0, NULL, NULL);

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

	// the extra data necessary for convolution
	cl_uint extra = ceil(FILTER_LENGTH/2);

	ncount_t * counts = createSampleCounts(NUM_NUCL + extra);
	printCounts(counts, 5);
	float * entropy_out = new float[NUM_NUCL + extra];
	float * convolution_out = new float[NUM_NUCL];
	cl_uint * threshold_out = new cl_uint[NUM_NUCL];
	
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

	/////////////////////////////// ENTROPY /////////////////////////////////
	counts_dev = clCreateBuffer(context, CL_MEM_READ_ONLY, NUM_NUCL * sizeof(ncount_t), NULL, &err);
	checkError(err);
	entropy_out_dev = clCreateBuffer(context, CL_MEM_WRITE_ONLY, NUM_NUCL * sizeof(float), NULL, &err);
	checkError(err);
	//printf("%d", sizeof(reference));
	clEnqueueWriteBuffer(command_queue, counts_dev, CL_TRUE, 0, NUM_NUCL * sizeof(ncount_t), counts, 0, NULL, NULL);

	program = build_program(context, device_id, num_of_devices, ENTROPY);

	kernel = clCreateKernel(program, "calculateEntropy", &err);
	checkError(err);
	for (int i = 1; i <= 1024; i = i*2) {
		cl_int grain_size = i;
		global = NUM_NUCL / GRAIN_SIZE;
		local = LOCAL;
		auto start_time = std::chrono::high_resolution_clock::now();
		runEntropyKernel(context, command_queue, program, num_of_devices, device_id, kernel, global, local, counts, entropy_out, grain_size, counts_dev, entropy_out_dev);
		auto current_time = std::chrono::high_resolution_clock::now();
		std::cout << "Grain Size: " << grain_size << " Time " << std::chrono::duration_cast<std::chrono::duration<float>>(current_time - start_time).count() << " seconds" << std::endl;
	}
	system("pause");
	/////////////////////////// CONVOLUTION /////////////////////////////
	float filter[] = { 1, 1, 1, 1, 1 };
	cl_uint filter_length = FILTER_LENGTH;
	cl_uint input_length = NUM_NUCL;
	global = NUM_NUCL;
	local = LOCAL;

	runConvolutionKernel(context, command_queue, program2, num_of_devices, device_id, kernel, global, local, filter, convolution_out, filter_length, input_length, filter_dev, convolution_out_dev, entropy_out_dev);
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

	runThresholdKernel(context, command_queue, program3, num_of_devices, device_id, kernel, global, local, threshold_out, threshold, threshold_out_dev, convolution_out_dev);
	
	delete []counts;
	delete []entropy_out;
	delete []convolution_out;
	delete []threshold_out;

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