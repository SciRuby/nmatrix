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
# SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
# NMatrix is Copyright (c) 2012, Ruby Science Foundation
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
# == nmatrix_spec.rb
#
# Element-wise operation tests.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")

describe NMatrix do

		it "should perform complex matrix operations" do
			x = Complex(3,4)
			y = Complex(1,2)			
			@n = NMatrix.new(:dense, 2, [x, 1, 1, y.conjugate], :complex64)
			@m = NMatrix.new(:dense, 2, [y, 1, 1, x.conjugate], :complex64)
			t = @n.dot(@m)
			t[0,0].should == Complex(-4,10)
			t[0,1].should == 6
			t[1,0].should == 2
			t[1,1].should == Complex(-4,-10)				
		end		
end
