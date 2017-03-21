cmake_minimum_required(VERSION 3.5)
project(kmers_2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    contig_generation.h
    kmer_hash.h
    packingDNAseq.h
    pgen_max_upc.cpp
    serial.c
    sort.cpp)

add_executable(kmers_2 ${SOURCE_FILES})