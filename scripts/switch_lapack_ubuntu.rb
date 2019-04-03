#!/usr/bin/env ruby

if ARGV[0] == "atlas"
  lapack_prefix = "/usr/lib/atlas-base/atlas"
  blas_prefix = "/usr/lib/atlas-base/atlas"
elsif ARGV[0] == "openblas"
  lapack_prefix = "/usr/lib/openblas-base"
  blas_prefix = "/usr/lib/openblas-base"
elsif ARGV[0] == "ref"
  lapack_prefix = "/usr/lib/lapack"
  blas_prefix = "/usr/lib/libblas"
else
  puts "options are atlas, openblas, or ref"
  exit
end

def run(cmd)
  puts "> #{cmd}"
  system cmd
end

run "update-alternatives --set liblapack.so.3 #{lapack_prefix}/liblapack.so.3"
run "update-alternatives --set liblapack.so #{lapack_prefix}/liblapack.so"
run "update-alternatives --set libblas.so.3 #{blas_prefix}/libblas.so.3"
run "update-alternatives --set libblas.so #{blas_prefix}/libblas.so"
