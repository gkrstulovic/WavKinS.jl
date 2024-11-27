#!/bin/bash

julia --threads 1 Bench_number_of_threads.jl
julia --threads 2 Bench_number_of_threads.jl
julia --threads 4 Bench_number_of_threads.jl
julia --threads 8 Bench_number_of_threads.jl
# julia --threads 16 Bench_number_of_threads.jl
# julia --threads 32 Bench_number_of_threads.jl
# julia --threads 64 Bench_number_of_threads.jl
# julia --threads 128 Bench_number_of_threads.jl