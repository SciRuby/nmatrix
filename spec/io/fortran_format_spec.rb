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
# Basic tests for NMatrix::IO::FortranFormat.
#

require './lib/nmatrix'

describe NMatrix::IO::FortranFormat do
  it "parses integer FORTRAN formats" do
    int_fmt =  NMatrix::IO::FortranFormat::Reader.new('(16I5)').parse

    expect(int_fmt[:format_code]).to eq "INT_ID"
    expect(int_fmt[:repeat])     .to eq 16
    expect(int_fmt[:field_width]).to eq 5    

    int_fmt = NMatrix::IO::FortranFormat::Reader.new('(I4)').parse 

    expect(int_fmt[:format_code]).to eq "INT_ID"
    expect(int_fmt[:field_width]).to eq 4
  end

  it "parses floating point FORTRAN formats" do
    fp_fmt = NMatrix::IO::FortranFormat::Reader.new('(10F7.1)').parse

    expect(fp_fmt[:format_code])       .to eq "FP_ID"
    expect(fp_fmt[:repeat])            .to eq 10
    expect(fp_fmt[:field_width])       .to eq 7
    expect(fp_fmt[:post_decimal_width]).to eq 1

    fp_fmt = NMatrix::IO::FortranFormat::Reader.new('(F4.2)').parse

    expect(fp_fmt[:format_code])       .to eq "FP_ID"
    expect(fp_fmt[:field_width])       .to eq 4
    expect(fp_fmt[:post_decimal_width]).to eq 2
  end

  it "parses exponential FORTRAN formats" do
	  exp_fmt = NMatrix::IO::FortranFormat::Reader.new('(2E8.3E3)').parse

	  expect(exp_fmt[:format_code])       .to eq "EXP_ID"
	  expect(exp_fmt[:repeat])            .to eq 2
	  expect(exp_fmt[:field_width])       .to eq 8
	  expect(exp_fmt[:post_decimal_width]).to eq 3
	  expect(exp_fmt[:exponent_width])    .to eq 3

	  exp_fmt = NMatrix::IO::FortranFormat::Reader.new('(3E3.6)').parse

	  expect(exp_fmt[:format_code])       .to eq "EXP_ID"
	  expect(exp_fmt[:repeat])            .to eq 3
	  expect(exp_fmt[:field_width])       .to eq 3
	  expect(exp_fmt[:post_decimal_width]).to eq 6

	  exp_fmt = NMatrix::IO::FortranFormat::Reader.new('(E4.5)').parse
	  expect(exp_fmt[:format_code])				.to eq "EXP_ID"
	  expect(exp_fmt[:field_width])				.to eq 4
	  expect(exp_fmt[:post_decimal_width]).to eq 5
  end

  ['I3', '(F4)', '(E3.', '(E4.E5)'].each do |bad_format|
	  it "doesn't let bad input through : #{bad_format}" do
		  expect {
			  NMatrix::IO::FortranFormat::Reader.new(bad_format).parse
		  }.to raise_error(IOError)
  end
end
end
