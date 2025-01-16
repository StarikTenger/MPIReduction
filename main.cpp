#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define VECTOR_SIZE 1000000000  // Size of the vector

int main(int argc, char *argv[]) {
    int rank, n_nodes;
    float *vector = NULL;
    float local_sum = 0.0, global_sum = 0.0;
    double start_time, end_time;

    MPI_Init(&argc, &argv);                // Initialize MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get process rank
    MPI_Comm_size(MPI_COMM_WORLD, &n_nodes); // Get total number of processes

    int chunk_size = VECTOR_SIZE / n_nodes;  // Divide vector among processes

    // Allocate memory for the local chunk of the vector
    float *local_vector = (float *)malloc(chunk_size * sizeof(float));

    // Initialize the vector on the root process
    if (rank == 0) {

        std::cout << "Number of processes: " << n_nodes << std::endl;
        vector = (float *)malloc(VECTOR_SIZE * sizeof(float));
        for (int i = 0; i < VECTOR_SIZE; i++) {
            vector[i] = 1.0;  // Fill the vector with 1.0 for simplicity
        }
        start_time = MPI_Wtime();  // Start the timer on the root process
    }

    // Scatter the vector to all processes
    MPI_Scatter(vector, chunk_size, MPI_FLOAT, local_vector, chunk_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Print the number of threads available for OpenMP
    #pragma omp parallel
    {
        #pragma omp master
        {
            int n_threads = omp_get_num_threads();
            int max_threads = omp_get_max_threads();
            char hostname[257];
            gethostname(hostname, sizeof hostname);
            printf("name=%s;\trank=%d;\topenMP threads: %d/%d\n", hostname, rank, n_threads, max_threads);
        }
    }

    // Parallel reduction within each process using OpenMP
    #pragma omp parallel for reduction(+:local_sum)
    for (int i = 0; i < chunk_size; i++) {
        local_sum += local_vector[i];
    }

    // Reduce local sums to the global sum using MPI
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Print the global sum and execution time on the root process
    if (rank == 0) {
        end_time = MPI_Wtime();  // Stop the timer on the root process
        printf("Global sum: %f\n", global_sum);
        printf("Execution time: %f seconds\n", end_time - start_time);
        free(vector);
    }

    free(local_vector);
    MPI_Finalize(); // Finalize MPI
    return 0;
}
