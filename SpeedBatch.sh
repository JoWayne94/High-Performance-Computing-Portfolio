#!/bin/sh

> executionTimes.txt
echo "Shell script for scaling plot" > executionTimes.txt

# 25 processes
echo "25 processes" >> executionTimes.txt
start=$(date +%s.%N)
mpiexec -np 25 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=5 --Py=5 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 25 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=5 --Py=5 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 25 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=5 --Py=5 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 25 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=5 --Py=5 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 25 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=5 --Py=5 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

# 16 processes
echo "16 processes" >> executionTimes.txt
start=$(date +%s.%N)
mpiexec -np 16 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=4 --Py=4 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 16 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=4 --Py=4 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 16 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=4 --Py=4 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 16 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=4 --Py=4 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 16 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=4 --Py=4 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

# 9 processes
echo "9 processes" >> executionTimes.txt
start=$(date +%s.%N)
mpiexec -np 9 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=3 --Py=3 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 9 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=3 --Py=3 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 9 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=3 --Py=3 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 9 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=3 --Py=3 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 9 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=3 --Py=3 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

# 4 processes
echo "4 processes" >> executionTimes.txt
start=$(date +%s.%N)
mpiexec -np 4 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=2 --Py=2 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 4 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=2 --Py=2 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 4 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=2 --Py=2 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 4 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=2 --Py=2 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 4 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=2 --Py=2 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

# 1 process
echo "1 process" >> executionTimes.txt
start=$(date +%s.%N)
mpiexec -np 1 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=1 --Py=1 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 1 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=1 --Py=1 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 1 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=1 --Py=1 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 1 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=1 --Py=1 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt

start=$(date +%s.%N)
mpiexec -np 1 ./Solve --Lx=1 --Ly=1 --Nx=161 --Ny=161 --Px=1 --Py=1 --dt=0.0005 --T=0.2 --Re=100
duration=$(echo "$(date +%s.%N) - $start" | bc)
execution_time=`printf "%.2f seconds" $duration`
echo "$execution_time" >> executionTimes.txt
