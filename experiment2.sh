#!/bin/bash 

#OAR -l nodes=2
#OAR -O OAR_%jobid%.out
#OAR -E OAR_%jobid%.err 

# display some information about attributed resources
hostname 
nvidia-smi 

echo "$(uniq $OAR_FILE_NODES | sed -n '1p') slots=1" > hostfile
echo "$(uniq $OAR_FILE_NODES | sed -n '2p') slots=1" >> hostfile

for floats in 1000 1000000 200000000; do
    filename="experiment2_${floats}_floats.txt"

    echo "Experiment 2: reduction on 2 nodes MPI + openMP (redefined reduce OP)" > $filename
    echo "" >> $filename

    for threads in 1 2 4 8 16 32 64; do
        echo " -------- Reduction with $threads threads --------" >> $filename
        echo "Experiment with $threads threads"
        mpirun -np 2 -hostfile hostfile -mca mtl psm2 -mca pml ^ucx,ofi -mca btl ^ofi,openib -x OMP_NUM_THREADS=$threads ~/MPIReduction/build/simple_reduce $floats >> $filename
        echo "" >> $filename
    done
done