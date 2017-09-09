//Includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <map>
#include <math.h>
#include "reads.h"
#include <chrono>
#include "mapping.h"
#include "C:\Users\korob\Desktop\common\basic.hpp"
#include "C:\Users\korob\Desktop\common\cmdparser.hpp"
#include "C:\Users\korob\Desktop\common\oclobject.hpp"


using namespace std;

#define KMER "kernel1.txt"
#define BSORT "bsort.txt"
#define FREQUENCY "frequency.txt"
#define REDUCTION "reduce_original.txt"
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <CL/cl.h>

#define LOCAL_SIZE 512 //512 independent threads per group
#define DATA_SIZE 4*LOCAL_SIZE //one work group
/* There are 64 * 512 reads in this sequence (1152 * 1024 bytes). Each thread is responsible for 64 reads*/
using namespace std;

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
		exit(1);
	}

	return program;
}

// OpenCL sort implementation
float ExecuteSortKernel(OpenCLBasic &oclobjects, OpenCLProgramOneKernel &executable, cl_uint* p_input, cl_int arraySize, cl_uint sortAscending)
{
	double   perf_start;
	double   perf_stop;

	cl_int err = CL_SUCCESS;
	cl_int numStages = 0;
	cl_uint temp;

	cl_int stage;
	cl_int passOfStage;

	// create OpenCL buffer using input array memory
	cl_mem cl_input_buffer = 
		clCreateBuffer
		(
			oclobjects.context,
			CL_MEM_USE_HOST_PTR,
			zeroCopySizeAlignment(sizeof(cl_int) * arraySize, oclobjects.device),
			p_input,
			&err
			);
	SAMPLE_CHECK_ERRORS(err);

	if (cl_input_buffer == (cl_mem)0)
	{
		throw Error("Failed to create input data Buffer\n");
	}

	for (temp = arraySize; temp > 2; temp >>= 1)
	{
		numStages++;
	}

	err = clSetKernelArg(executable.kernel, 0, sizeof(cl_mem), (void *)&cl_input_buffer);
	SAMPLE_CHECK_ERRORS(err);
	err = clSetKernelArg(executable.kernel, 3, sizeof(cl_uint), (void *)&sortAscending);
	SAMPLE_CHECK_ERRORS(err);

	perf_start = time_stamp();
	for (stage = 0; stage < numStages; stage++)
	{
		err = clSetKernelArg(executable.kernel, 1, sizeof(cl_uint), (void *)&stage);
		SAMPLE_CHECK_ERRORS(err);

		for (passOfStage = stage; passOfStage >= 0; passOfStage--)
		{
			err = clSetKernelArg(executable.kernel, 2, sizeof(cl_uint), (void *)&passOfStage);
			SAMPLE_CHECK_ERRORS(err);

			// set work-item dimensions
			size_t gsz = arraySize / (2 * 4);
			size_t global_work_size[1] = { passOfStage ? gsz : gsz << 1 };    //number of quad items in input array

																			  // execute kernel
			err = clEnqueueNDRangeKernel(oclobjects.queue, executable.kernel, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
			SAMPLE_CHECK_ERRORS(err);
		}
	}

	err = clFinish(oclobjects.queue);
	SAMPLE_CHECK_ERRORS(err);

	perf_stop = time_stamp();

	void* tmp_ptr = NULL;
	tmp_ptr = clEnqueueMapBuffer(oclobjects.queue, cl_input_buffer, true, CL_MAP_READ, 0, sizeof(cl_int) * arraySize, 0, NULL, NULL, &err);
	SAMPLE_CHECK_ERRORS(err);
	if (tmp_ptr != p_input)
	{
		throw Error("clEnqueueMapBuffer failed to return original pointer\n");
	}

	err = clFinish(oclobjects.queue);
	SAMPLE_CHECK_ERRORS(err);

	err = clEnqueueUnmapMemObject(oclobjects.queue, cl_input_buffer, tmp_ptr, 0, NULL, NULL);
	SAMPLE_CHECK_ERRORS(err);

	err = clReleaseMemObject(cl_input_buffer);
	SAMPLE_CHECK_ERRORS(err);

	return (float)(perf_stop - perf_start);
}

int main(int argc, const char** argv)
{
	CmdParserCommon cmdparser(argc, argv);

	cl_context context;
	cl_context_properties properties[3];
	cl_kernel kernel;
	cl_command_queue command_queue;
	cl_program program, program1, program2, program3;
	cl_int err;
	cl_uint num_of_platforms = 0;
	cl_platform_id platform_id;
	cl_device_id device_id;
	cl_uint num_of_devices = 0;
	cl_mem reads, output, sorted_output, in, out, out1, final_in, final, freq_index;

	size_t global, local;

	size_t numWorkGroups = DATA_SIZE / LOCAL_SIZE;
	int k = 8;
	unsigned int kmers = (128 - k + 1) * 32 * 1024;
	unsigned int num_kmers = nextPowerOf2(kmers);
	size_t resultSize = num_kmers*4;
	cout << "total elements: " << num_kmers << endl << "kmers: " << kmers << endl;
	unsigned int * sample = createSampleReads(LENGTH);
	cl_uint *table = new unsigned int[num_kmers];
	int sequence_length = 1152 * 1024;

	cl_queue_properties props[] = {
		CL_QUEUE_PROPERTIES,
		CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_ON_DEVICE | CL_QUEUE_ON_DEVICE_DEFAULT,
		0
	};

	// retreive a list of platforms avaible
	if (clGetPlatformIDs(1, &platform_id, &num_of_platforms) != CL_SUCCESS)
	{
		printf("Unable to get platform_id\n");
		return 1;
	}

	// try to get a supported GPU device
	if (clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &num_of_devices) != CL_SUCCESS)
	{
		printf("Unable to get device_id\n");
		return 1;
	}

	// context properties list - must be terminated with 0
	properties[0] = CL_CONTEXT_PLATFORM;
	properties[1] = (cl_context_properties)platform_id;
	properties[2] = 0;

	// create a context with the GPU device
	context = clCreateContext(properties, 1, &device_id, NULL, NULL, &err);

	// create command queue using the context and device
	command_queue = clCreateCommandQueueWithProperties(context, device_id, props, &err);


	program = build_program(context, device_id, num_of_devices, KMER);

	reads = clCreateBuffer(context, CL_MEM_READ_ONLY, sequence_length, NULL, &err);
	output = clCreateBuffer(context, CL_MEM_WRITE_ONLY, resultSize, NULL, &err);

	// load data into the input buffer
	clEnqueueWriteBuffer(command_queue, reads, CL_TRUE, 0, sequence_length, sample, 0, NULL, NULL);
	//clEnqueueWriteBuffer(command_queue, output, CL_TRUE, 0, resultSize, result, 0, NULL, NULL);

	// specify which kernel from the program to execute
	kernel = clCreateKernel(program, "createMap", &err);
	//printf("%d\n", sizeof(cl_mem));
	// set the argument list for the kernel command
	err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &reads);
	err |= clSetKernelArg(kernel, 1, sizeof(unsigned int), &k);
	err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &output);
	err |= clSetKernelArg(kernel, 3, sizeof(unsigned int), &sequence_length);

	if (err != CL_SUCCESS) {
		printf("Unable to set args");
		return 1;
	}
	global = DATA_SIZE;
	local = LOCAL_SIZE;
	auto start_time = std::chrono::high_resolution_clock::now();
	// enqueue the kernel command for execution
	err = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global, &local, 0, NULL, NULL);
	if (err != CL_SUCCESS) {
		printf("something wrong %d", err);
		system("pause");
		return 1;
	}
	clFinish(command_queue);

	auto current_time = std::chrono::high_resolution_clock::now();
	std::cout << "Program has been running for " << std::chrono::duration_cast<std::chrono::duration<float>>(current_time - start_time).count() << " seconds" << std::endl;
	// copy the results from out of the output buffer
	clEnqueueReadBuffer(command_queue, output, CL_TRUE, 0, resultSize, table, 0, NULL, NULL);

	print_elements(&table[0], 30);
	//system("pause");
	// bitonic sort /////////////////////////////////////////////////////////////////////////////////////////////////

	CmdOption<bool> reverse_sort(cmdparser, 0, "reverse-sort", "", "performs descending sort (default is ascending).", false);
	CmdOption<int> array_size(cmdparser, 's', "size", "<integer>", "input/output array size.", num_kmers);
	cmdparser.parse();
	// Create the necessary OpenCL objects up to device queue.
	OpenCLBasic oclobjects(
		cmdparser.platform.getValue(),
		cmdparser.device_type.getValue(),
		cmdparser.device.getValue()
		);
	printf("Sort order is %s\n", reverse_sort.getValue() ? "descending" : "ascending");
	printf("Input size is %d items\n", array_size.getValue());

	cl_uint dev_alignment = zeroCopyPtrAlignment(oclobjects.device);
	size_t aligned_size = zeroCopySizeAlignment(sizeof(cl_int) * array_size.getValue(), oclobjects.device);
	//printf("OpenCL data alignment is %d bytes.\n", dev_alignment);
	//cl_int* p_input = (cl_int*)aligned_malloc(aligned_size, dev_alignment);
	//generateInput(p_input, array_size.getValue());
	// Build kernel
	OpenCLProgramOneKernel executable(oclobjects, L"bitonicsort.txt", "", "BitonicSort");

	float ocl_time = ExecuteSortKernel(oclobjects, executable, table, array_size.getValue(), !reverse_sort.getValue());
	print_elements(&table[num_kmers-kmers] , 300);
	//system("pause");
	///////////// frequency ////////////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int f_index = 0;
	kmer_t *frequency_table = new kmer_t[num_kmers];
	kmer_t *frequency_table1 = new kmer_t[num_kmers];
	kmer_t *frequency_table_final = new kmer_t[num_kmers];
	intptr_t *tables_in = new intptr_t[100];
	program2 = build_program(context, device_id, num_of_devices, FREQUENCY);

	in = clCreateBuffer(context, CL_MEM_READ_ONLY, num_kmers*4, NULL, &err);
	final_in = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(tables_in)*4, NULL, &err);
	out = clCreateBuffer(context, CL_MEM_READ_WRITE, num_kmers*sizeof(kmer_t), NULL, &err);
	out1 = clCreateBuffer(context, CL_MEM_READ_WRITE, num_kmers*sizeof(kmer_t), NULL, &err);
	freq_index = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(unsigned int), NULL, &err);

	clEnqueueWriteBuffer(command_queue, freq_index, CL_TRUE, 0, sizeof(unsigned int), &f_index, 0, NULL, NULL);
	clEnqueueWriteBuffer(command_queue, in, CL_TRUE, 0, num_kmers*4, table, 0, NULL, NULL);

	kernel = clCreateKernel(program2, "countFrequency", &err);

	// set the argument list for the kernel command
	err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &in);
	err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &out);
	err |= clSetKernelArg(kernel, 2, sizeof(unsigned int), &num_kmers);
	err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &final_in);
	err |= clSetKernelArg(kernel, 4, sizeof(unsigned int), &freq_index);
	if (err != CL_SUCCESS) {
		printf("Unable to set args");
		system("pause");
		return 1;
	}
	// enqueue the kernel command for execution...second table
	local = 1;
	global = 1;
	err = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global, &local, 0, NULL, NULL);
	if (err != CL_SUCCESS) {
		printf("something wrong %d", err);
		system("pause");
		return 1;
	}
	clFinish(command_queue);
	// set the argument list for the kernel command
	err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &in);
	err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &out1);
	err |= clSetKernelArg(kernel, 2, sizeof(unsigned int), &num_kmers);
	err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &final_in);
	err |= clSetKernelArg(kernel, 4, sizeof(unsigned int), &freq_index);
	if (err != CL_SUCCESS) {
		printf("Unable to set args");
		system("pause");
		return 1;
	}

	err = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global, &local, 0, NULL, NULL);
	if (err != CL_SUCCESS) {
		printf("something wrong %d", err);
		system("pause");
		return 1;
	}
	clFinish(command_queue);

	//std::cout << "Program has been running for " << std::chrono::duration_cast<std::chrono::duration<float>>(current_time - start_time).count() << " seconds" << std::endl;
	// copy the results from out of the output buffer
	clEnqueueReadBuffer(command_queue, out, CL_TRUE, 0, num_kmers*sizeof(kmer_t), frequency_table, 0, NULL, NULL);
	clEnqueueReadBuffer(command_queue, out1, CL_TRUE, 0, num_kmers*sizeof(kmer_t), frequency_table1, 0, NULL, NULL);
	print_freq_table(frequency_table, 10);
	//system("pause");
	print_freq_table(frequency_table1, 10);
	system("pause");
	///////////////// reduction ///////////////////////////////////////////////////////////////////////////////////////
	program3 = build_program(context, device_id, num_of_devices, REDUCTION);
	int num_tables = 2;
	local = 2;
	global = 2;
	unsigned int flag = 0;
	unsigned int out_ind = 0;

	cout << sizeof(tables_in) << endl;
	final = clCreateBuffer(context, CL_MEM_READ_WRITE, num_kmers*sizeof(kmer_t), NULL, &err);
	kernel = clCreateKernel(program3, "reduceFrequency", &err);

	// set the argument list for the kernel command
	err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &final_in);
	if (err != CL_SUCCESS) {
		printf("1. Unable to set args %d\n", err);
		system("pause");
		return 1;
	}
	err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &final);
	if (err != CL_SUCCESS) {
		printf("2. Unable to set args %d\n", err);
		system("pause");
		return 1;
	}
	err |= clSetKernelArg(kernel, 2, sizeof(unsigned int), &num_kmers);
	err |= clSetKernelArg(kernel, 3, sizeof(unsigned int), &out_ind);
	err |= clSetKernelArg(kernel, 4, sizeof(unsigned int), &flag);
	if (err != CL_SUCCESS) {
		printf("Unable to set args %d\n", err);
		system("pause");
		return 1;
	}
	// enqueue the kernel command for execution
	err = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global, &local, 0, NULL, NULL);
	if (err != CL_SUCCESS) {
		printf("something wrong %d\n", err);
		system("pause");
		return 1;
	}
	clFinish(command_queue); 
	clEnqueueReadBuffer(command_queue, final, CL_TRUE, 0, num_kmers*sizeof(kmer_t), frequency_table_final, 0, NULL, NULL);
	print_freq_table(frequency_table_final, 10);
	system("pause");
	// cleanup - release OpenCL resources
	clReleaseMemObject(reads);
	clReleaseMemObject(output);
	clReleaseProgram(program);
	clReleaseKernel(kernel);
	clReleaseCommandQueue(command_queue);
	clReleaseContext(context);

	delete[] sample;
	delete[] table;
	delete[] frequency_table;
	return 0;
}