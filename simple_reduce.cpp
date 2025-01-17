#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <random>

void vector_sum(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
    float *in = (float *)invec;
    float *inout = (float *)inoutvec;

    #pragma omp parallel for
    for (int i = 0; i < *len; i++) {
        inout[i] += in[i];
    }
}

int main(int argc, char *argv[]) {
    int rank, n_nodes;
    double start_time, scatter_time, end_time;

    MPI_Init(&argc, &argv);                // Initialize MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get process rank
    MPI_Comm_size(MPI_COMM_WORLD, &n_nodes); // Get total number of processes

    size_t VECTOR_SIZE;
    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <vector_size>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }
    VECTOR_SIZE = atoi(argv[1]);
    
    float *local_vector = new float[VECTOR_SIZE];

    #pragma omp parallel
    {
        std::random_device rd;
        std::mt19937 gen(rd() + omp_get_thread_num() + rank * 1000);
        std::uniform_real_distribution<float> dis(-100.0, 100.0);
        
        #pragma omp for
        for (int i = 0; i < VECTOR_SIZE; i++) {
            local_vector[i] = dis(gen);
        }
    }

    if (rank == 0) {
        std::cout << "Number of processes: " << n_nodes << std::endl;
        std::cout << "Vector size: " << VECTOR_SIZE << std::endl;
    }

    // Print the number of threads available for OpenMP
    #pragma omp parallel
    {
        #pragma omp master
        {
            int n_threads = omp_get_num_threads();
            int max_threads = omp_get_max_threads();
            char hostname[257];
            gethostname(hostname, sizeof(hostname));
            std::cout << "name=" << hostname << ";\trank=" << rank << ";\topenMP threads: " << n_threads << "/" << max_threads << std::endl;
        }
    }

    float *vector = new float[VECTOR_SIZE];

    MPI_Op vector_sum_op;
    MPI_Op_create(&vector_sum, 1, &vector_sum_op);

    start_time = MPI_Wtime();

    MPI_Reduce(local_vector, vector, VECTOR_SIZE, MPI_FLOAT, vector_sum_op, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        end_time = MPI_Wtime();
        std::cout << "Elapsed time: " << (end_time - start_time) * 1000 << " ms" << std::endl;
    }

    delete[] vector;
    delete[] local_vector;

    MPI_Finalize(); // Finalize MPI
}
