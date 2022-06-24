cmake_minimum_required(VERSION 3.10)

execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

project(vgrna-project-scripts)

ExternalProject_Add(seqlib-proj
  GIT_REPOSITORY "https://github.com/walaj/SeqLib.git"
  GIT_TAG "08771285f60f3434c796cd4126ea02c577632ddb"
  SOURCE_DIR "${PROJECT_SOURCE_DIR}/external/seqlib"
  BUILD_IN_SOURCE TRUE
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ./configure
  INSTALL_COMMAND ""
)

include_directories(
  ${PROJECT_SOURCE_DIR}/external/seqlib
  ${PROJECT_SOURCE_DIR}/external/seqlib/htslib
)

link_directories(
  ${PROJECT_SOURCE_DIR}/external/seqlib/src
  ${PROJECT_SOURCE_DIR}/external/seqlib/htslib
  ${PROJECT_SOURCE_DIR}/external/seqlib/bwa
)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  set(CMAKE_CXX_FLAGS "--std=c++14 -Xpreprocessor -fopenmp -g -O3 -Wall -DGIT_COMMIT='\"${GIT_COMMIT}\"'")

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

  set(CMAKE_CXX_FLAGS "--std=c++14 -fopenmp -g -O3 -DGIT_COMMIT='\"${GIT_COMMIT}\"'")

endif()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_executable(calc_allele_read_coverage src/calc_allele_read_coverage.cpp)
add_dependencies(calc_allele_read_coverage seqlib-proj)
target_link_libraries(calc_allele_read_coverage seqlib hts z)

add_executable(calc_exon_read_coverage src/calc_exon_read_coverage.cpp)
add_dependencies(calc_exon_read_coverage seqlib-proj)
target_link_libraries(calc_exon_read_coverage seqlib hts z)

add_executable(convert_rsem_sim_info src/convert_rsem_sim_info.cpp)
add_dependencies(convert_rsem_sim_info seqlib-proj)
target_link_libraries(convert_rsem_sim_info seqlib hts z)

add_executable(calc_vg_benchmark_stats src/calc_vg_benchmark_stats.cpp)
add_dependencies(calc_vg_benchmark_stats seqlib-proj)
target_link_libraries(calc_vg_benchmark_stats seqlib hts z)

add_executable(calc_read_regions_overlap_stats src/calc_read_regions_overlap_stats.cpp)
add_dependencies(calc_read_regions_overlap_stats seqlib-proj)
target_link_libraries(calc_read_regions_overlap_stats seqlib hts z)

add_executable(calc_read_transcript_overlap_stats src/calc_read_transcript_overlap_stats.cpp)
add_dependencies(calc_read_transcript_overlap_stats seqlib-proj)
target_link_libraries(calc_read_transcript_overlap_stats seqlib hts z)

add_executable(convert_snaptron_to_bed src/convert_snaptron_to_bed.cpp)
add_dependencies(convert_snaptron_to_bed seqlib-proj)
target_link_libraries(convert_snaptron_to_bed seqlib hts z)

add_executable(filter_multi_gene_introns src/filter_multi_gene_introns.cpp)
add_dependencies(filter_multi_gene_introns seqlib-proj)
target_link_libraries(filter_multi_gene_introns seqlib hts z)

add_executable(convert_cds_alignments_to_haplotypes src/convert_cds_alignments_to_haplotypes.cpp)
add_dependencies(convert_cds_alignments_to_haplotypes seqlib-proj)
target_link_libraries(convert_cds_alignments_to_haplotypes seqlib hts bwa z)
