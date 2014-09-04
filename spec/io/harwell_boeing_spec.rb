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
# == io_spec.rb
#
# Basic tests for NMatrix::IO::HarwelBoeing.

# TODO : After the fortran format thing is done
require 'spec_helper'
require "./lib/nmatrix"

describe NMatrix::IO::HarwellBoeing do
  def check_file_header header
    expect(header[:title])    .to eq("Title")
    expect(header[:key])      .to eq("Key")

    expect(header[:totcrd])   .to eq(5)
    expect(header[:ptrcrd])   .to eq(1)
    expect(header[:indcrd])   .to eq(1)
    expect(header[:valcrd])   .to eq(3)
    expect(header[:rhscrd])   .to eq(0)
    
    expect(header[:mxtype])   .to eq('RUA')
    expect(header[:nrow])     .to eq(5)
    expect(header[:ncol])     .to eq(5)
    expect(header[:nnzero])   .to eq(13)
    expect(header[:neltvl])   .to eq(0)

    expect(header[:ptrfmt])   .to eq({
      format_code: "INT_ID",
      repeat:             6,         
      field_width:        3
      })
    expect(header[:indfmt])   .to eq({
      format_code: "INT_ID",
      repeat:            13,
      field_width:       3
      })
    expect(header[:valfmt])   .to eq({
      format_code:         "EXP_ID",
      repeat:                     5,
      field_width:               15,
      post_decimal_width:         8
      })
    expect(header[:rhsfmt])   .to eq({
      format_code:         "EXP_ID",
      repeat:                     5,
      field_width:               15,
      post_decimal_width:         8
      })
  end

  it "loads a Harwell Boeing file values and header (currently real only)" do
    n, h = NMatrix::IO::HarwellBoeing.load("spec/io/test.rua")

    expect(n.is_a? NMatrix).to eq(true)
    expect(n.cols)         .to eq(5)
    expect(n.rows)         .to eq(5)

    expect(n[0,0])         .to eq(11)
    expect(n[4,4])         .to eq(55)

    expect(h.is_a? Hash).to eq(true) 
    check_file_header(h)
  end

  it "loads only the header of the file when specified" do
    h = NMatrix::IO::HarwellBoeing.load("spec/io/test.rua", header: true)

    expect(h.is_a? Hash).to eq(true)
    check_file_header(h)
  end

  it "raises error for wrong Harwell Boeing file name" do
    expect{
      NMatrix::IO::HarwellBoeing.load("spec/io/wrong.afx")
    }.to raise_error(IOError)
  end
end