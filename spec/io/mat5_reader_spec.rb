# = NMatrix
#
# A linear algebra library for scientific computation in Ruby.
# NMatrix is part of SciRuby.
#
# NMatrix was originally inspired by and derived from NArray, by
# Masahiro Tanaka: http://narray.rubyforge.org
#
# == Copyright Information
#
# SciRuby is Copyright (c) 2010 - 2014, Ruby Science Foundation
# NMatrix is Copyright (c) 2012 - 2014, John Woods and the Ruby Science Foundation
#
# Please see LICENSE.txt for additional copyright notices.
#
# == Contributing
#
# By contributing source code to SciRuby, you agree to be bound by
# our Contributor Agreement:
#
# * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
#
# == fortran_format_spec.rb
#
# Basic tests for NMatrix::IO::Matlab
#

require_relative '../../lib/nmatrix'

describe NMatrix::IO::Matlab do
  def check_matrix_data_sparse arr
    expect(arr[0][:nonzero_max]).to eq(11)
    expect(arr[0][:matlab_class]).to eq(:mxSPARSE)
    expect(arr[0][:dimensions]).to eq([4,4])
    expect(arr[0][:matlab_name]).to eq("x")
    expect(arr[0][:real_part][:tag][:data_type]).to eq(:miDOUBLE)
    expect(arr[0][:real_part][:tag][:raw_data_type]).to eq(9)
    expect(arr[0][:real_part][:tag][:bytes]).to eq(32)
    expect(arr[0][:row_index][:tag][:data_type]).to eq(:miINT32)
    expect(arr[0][:row_index][:tag][:raw_data_type]).to eq(5)
    expect(arr[0][:row_index][:tag][:bytes]).to eq(20)
    expect(arr[0][:column_index][:tag][:data_type]).to eq(:miINT32)
    expect(arr[0][:column_index][:tag][:raw_data_type]).to eq(5)
    expect(arr[0][:column_index][:tag][:bytes]).to eq(44)
  end
  def check_matrix_data_dense arr
    expect(arr[0][:nonzero_max]).to eq(0)
    expect(arr[0][:matlab_class]).to eq(:mxDOUBLE)
    expect(arr[0][:dimensions]).to eq([4,5])
    expect(arr[0][:matlab_name]).to eq("x")
    expect(arr[0][:real_part][:tag][:data_type]).to eq(:miUINT8)
    expect(arr[0][:real_part][:tag][:raw_data_type]).to eq(2)
    expect(arr[0][:real_part][:tag][:bytes]).to eq(20)
  end
  it "reads MATLAB .mat file containing a single square sparse matrix" do
    pending("not yet implemented for NMatrix-JRuby") if jruby?
    n = NMatrix::IO::Matlab.load_mat("spec/4x4_sparse.mat")
    expect(n[0,0]).to eq(2)
    expect(n[1,1]).to eq(3)
    expect(n[1,3]).to eq(5)
    expect(n[3,0]).to eq(4)
    expect(n[2,2]).to eq(0)
    expect(n[3,3]).to eq(0)
  end

  it "reads MATLAB .mat file containing a single dense integer matrix" do
    n = NMatrix::IO::Matlab.load_mat("spec/4x5_dense.mat")
    m = NMatrix.new([4,5], [16,17,18,19,20,15,14,13,12,11,6,7,8,9,10,5,4,3,2,1])
    expect(n).to eq(m)
  end

  it "reads MATLAB .mat file containing a single dense double matrix" do
    n = NMatrix::IO::Matlab.load_mat("spec/2x2_dense_double.mat")
    m = NMatrix.new(2, [1.1, 2.0, 3.0, 4.0], dtype: :float64)
    expect(n).to eq(m)
  end
  it "reads MATLAB .mat file containing a single square sparse matrix111" do
    m =  NMatrix::IO::Matlab::Mat5Reader.new(File.open("spec/4x4_sparse.mat", "rb+"))
    n = m.to_ruby
    expect(m.byte_order).to eq(:little)
    expect(m.file_header.version).to eq(256)
    expect(m.file_header.endian).to eq("IM")
    expect(m.file_header.desc).to eq("MATLAB 5.0 MAT-file, Platform: GLNX86, Created on: Wed Apr 18 13:43:38 2012")
    expect(n[0,0]).to eq(2)
    expect(n[1,1]).to eq(3)
    expect(n[1,3]).to eq(5)
    expect(n[2,2]).to eq(0)
    expect(n[3,3]).to eq(0)
    check_matrix_data_sparse(m.to_a)
  end
  it "reads MATLAB .mat file containing a single dense double matrix123" do
    m =  NMatrix::IO::Matlab::Mat5Reader.new(File.open("spec/4x5_dense.mat", "rb+"))
    n = m.to_ruby
    expect(m.byte_order).to eq(:little)
    expect(m.file_header.version).to eq(256)
    expect(m.file_header.endian).to eq("IM")
    expect(m.file_header.desc).to eq("MATLAB 5.0 MAT-file, Platform: GLNX86, Created on: Wed Apr 18 19:17:17 2012")
    p = NMatrix.new([4,5], [16,17,18,19,20,15,14,13,12,11,6,7,8,9,10,5,4,3,2,1])
    expect(n).to eq(p)
    check_matrix_data_dense(m.to_a)
  end

end
