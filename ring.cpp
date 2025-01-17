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

void my_reduce(float* local_vector, float* vector, size_t VECTOR_SIZE, size_t rank, size_t n_nodes) {
    float *recvbuf = new float[VECTOR_SIZE / n_nodes];

    int chunk_size = VECTOR_SIZE / n_nodes;

    for (int step = 1; step < n_nodes; step++) {
        int send_rank = (rank - 1 + n_nodes) % n_nodes;
        int recv_rank = (rank + 1) % n_nodes;

        MPI_Sendrecv(local_vector + chunk_size * ((rank + step) % n_nodes), chunk_size, MPI_FLOAT, send_rank, 0,
                 recvbuf, chunk_size, MPI_FLOAT, recv_rank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        #pragma omp parallel for
        for (int i = 0; i < chunk_size; i++) {
            local_vector[chunk_size * ((rank + step + 1) % n_nodes) + i] += recvbuf[i];
        }
    }


    MPI_Gather(local_vector + rank * chunk_size, chunk_size, MPI_FLOAT, vector, chunk_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    delete[] recvbuf;
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
            local_vector[i] = 1.;//dis(gen);
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

    // ------------ Canonical reduce algorithm ------------

    // MPI_Reduce(local_vector, vector, VECTOR_SIZE, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    // if (rank == 0) {
    //     std::cout << "Canonical vector sum[0]: " << vector[0] << std::endl;
    // }

    // -----------------------------------------------------------


    MPI_Op vector_sum_op;
    MPI_Op_create(&vector_sum, 1, &vector_sum_op);

    start_time = MPI_Wtime();

    my_reduce(local_vector, vector, VECTOR_SIZE, rank, n_nodes);

    if (rank == 0) {
        end_time = MPI_Wtime();
        std::cout << "Elapsed time: " << (end_time - start_time) * 1000 << " ms" << std::endl;
        
        // std::cout << "Vector sum: ";
        // for (size_t i = 0; i < VECTOR_SIZE; i++) {
        //     std::cout << vector[i] << " ";
        // }
        // std::cout << std::endl;
    }


    

    MPI_Op_free(&vector_sum_op);

    delete[] vector;
    delete[] local_vector;

    MPI_Finalize(); // Finalize MPI
}
