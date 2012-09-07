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
# Basic tests for NMatrix.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")

describe NMatrix do

  it "calculates exact determinants on small square matrices" do
    a = NMatrix.new(:dense, 2, [1,2,3,4], :int64)
    x = a.det_exact
  end

  it "allows stype casting of a dim 2 matrix between dense, sparse, and list (different dtypes)" do
    m = NMatrix.new(:dense, [3,3], [0,0,1,0,2,0,3,4,5], :int64).
        cast(:yale, :int32).
        cast(:dense, :float64).
        cast(:list, :int32).
        cast(:dense, :int16).
        cast(:list, :int32).
        cast(:yale, :int64) #.
        #cast(:list, :int32).
        #cast(:dense, :int16)
    #m.should.equal?(original)
    # For some reason this causes some weird garbage collector problems when we uncomment these. The above lines won't
    # work at all in IRB, but work fine when run in a regular Ruby session.
  end


  it "compares two list matrices" do
    n = NMatrix.new(:list, [3,3,3], :int64)
    m = NMatrix.new(:list, [3,3,3], :int64)
    n.should == m
    n[0,0,0] = 5
    n.should_not == m
    n[0,0,1] = 52
    n[1,2,1] = -4

    m[0,0,0] = 5
    m[0,0,1] = 52
    m[1,2,1] = -4
    n.should == m
  end

  it "fills dense Ruby object matrix with nil" do
    n = NMatrix.new([4,3], :object)
    n[0,0].should == nil
  end

  it "fills dense with individual assignments" do
    n = NMatrix.new([4,3], :float64)
    n[0,0] = 14.0
    n[0,1] = 9.0
    n[0,2] = 3.0
    n[1,0] = 2.0
    n[1,1] = 11.0
    n[1,2] = 15.0
    n[2,0] = 0.0
    n[2,1] = 12.0
    n[2,2] = 17.0
    n[3,0] = 5.0
    n[3,1] = 2.0
    n[3,2] = 3.0

    n[0,0].should == 14.0
    n[0,1].should == 9.0
    n[0,2].should == 3.0
    n[1,0].should == 2.0
    n[1,1].should == 11.0
    n[1,2].should == 15.0
    n[2,0].should == 0.0
    n[2,1].should == 12.0
    n[2,2].should == 17.0
    n[3,0].should == 5.0
    n[3,1].should == 2.0
    n[3,2].should == 3.0
  end

  it "fills dense with a single mass assignment" do
    n = NMatrix.new([4,3], [14.0, 9.0, 3.0, 2.0, 11.0, 15.0, 0.0, 12.0, 17.0, 5.0, 2.0, 3.0])

    n[0,0].should == 14.0
    n[0,1].should == 9.0
    n[0,2].should == 3.0
    n[1,0].should == 2.0
    n[1,1].should == 11.0
    n[1,2].should == 15.0
    n[2,0].should == 0.0
    n[2,1].should == 12.0
    n[2,2].should == 17.0
    n[3,0].should == 5.0
    n[3,1].should == 2.0
    n[3,2].should == 3.0
  end

  it "fills dense with a single mass assignment, with dtype specified" do
    m = NMatrix.new([4,3], [14.0, 9.0, 3.0, 2.0, 11.0, 15.0, 0.0, 12.0, 17.0, 5.0, 2.0, 3.0], :float32)
    m[0,0].should == 14.0
    m[0,1].should == 9.0
    m[0,2].should == 3.0
    m[1,0].should == 2.0
    m[1,1].should == 11.0
    m[1,2].should == 15.0
    m[2,0].should == 0.0
    m[2,1].should == 12.0
    m[2,2].should == 17.0
    m[3,0].should == 5.0
    m[3,1].should == 2.0
    m[3,2].should == 3.0
  end


  it "list handles missing initialization value" do
    NMatrix.new(:list, 3, :int8)[0,0].should    == 0
    NMatrix.new(:list, 4, :float64)[0,0].should == 0.0
  end


  ##TODO: Make this test better. It's not nearly exhaustive enough as is.
  it "list handles recursive removal" do
    n = NMatrix.new(:list, [3,3,3], 0)
    n[0,0,0] = 2
    n[1,1,1] = 1
    n[1,0,0] = 3
    n[0,0,1] = 4

    n[0,0,0].should == 2
    n[1,1,1].should == 1
    n[1,0,0].should == 3
    n[0,0,1].should == 4

    n[1,1,1] = 0
    n[0,0,0].should == 2
    n[1,1,1].should == 0
    n[1,0,0].should == 3
    n[0,0,1].should == 4
  end


  it "dense handles missing initialization value" do
    n = NMatrix.new(3, :int8)
    n.stype.should == :dense
    n.dtype.should == :int8

    m = NMatrix.new(4, :float64)
    m.stype.should == :dense
    m.dtype.should == :float64
  end

  it "dense pretty_prints complex values" do
    n = NMatrix.new([4,3], COMPLEX_MATRIX43A_ARRAY, :complex128)
    n.pretty_print
  end

   # TODO: Get it working with ROBJ too
  [:byte,:int8,:int16,:int32,:int64,:float32,:float64,:rational64,:rational128].each do |left_dtype|
    [:byte,:int8,:int16,:int32,:int64,:float32,:float64,:rational64,:rational128].each do |right_dtype|

      # Won't work if they're both 1-byte, due to overflow.
      next if [:byte,:int8].include?(left_dtype) && [:byte,:int8].include?(right_dtype)

      # For now, don't bother testing int-int mult.
      #next if [:int8,:int16,:int32,:int64].include?(left_dtype) && [:int8,:int16,:int32,:int64].include?(right_dtype)
      it "dense handles #{left_dtype.to_s} dot #{right_dtype.to_s} matrix multiplication" do
        #STDERR.puts "dtype=#{dtype.to_s}"
        #STDERR.puts "2"

        nary = if left_dtype.to_s =~ /complex/
                 COMPLEX_MATRIX43A_ARRAY
               elsif left_dtype.to_s =~ /rational/
                 RATIONAL_MATRIX43A_ARRAY
               else
                 MATRIX43A_ARRAY
               end

        mary = if right_dtype.to_s =~ /complex/
                 COMPLEX_MATRIX32A_ARRAY
               elsif right_dtype.to_s =~ /rational/
                 RATIONAL_MATRIX32A_ARRAY
               else
                 MATRIX32A_ARRAY
               end

        n = NMatrix.new([4,3], nary, left_dtype)
        m = NMatrix.new([3,2], mary, right_dtype)

        m.shape[0].should == 3
        m.shape[1].should == 2
        m.dim.should == 2

        n.shape[0].should == 4
        n.shape[1].should == 3
        n.dim.should == 2

        n.shape[1].should == m.shape[0]

        r = n.dot m

        r[0,0].should == 273.0
        r[0,1].should == 455.0
        r[1,0].should == 243.0
        r[1,1].should == 235.0
        r[2,0].should == 244.0
        r[2,1].should == 205.0
        r[3,0].should == 102.0
        r[3,1].should == 160.0

        #r.dtype.should == :float64 unless left_dtype == :float32 && right_dtype == :float32
      end
    end
  end

  [:byte,:int8,:int16,:int32,:int64,:float32,:float64,:rational64,:rational128].each do |left_dtype|
    [:byte,:int8,:int16,:int32,:int64,:float32,:float64,:rational64,:rational128].each do |right_dtype|

      # Won't work if they're both 1-byte, due to overflow.
      next if [:byte,:int8].include?(left_dtype) && [:byte,:int8].include?(right_dtype)

      it "dense handles #{left_dtype.to_s} dot #{right_dtype.to_s} vector multiplication" do
        #STDERR.puts "dtype=#{dtype.to_s}"
        #STDERR.puts "2"
        n = NMatrix.new([4,3], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0], left_dtype)

        m = NVector.new(3, [2.0, 1.0, 0.0], right_dtype)

        m.shape[0].should == 3
        m.shape[1].should == 1

        n.shape[0].should == 4
        n.shape[1].should == 3
        n.dim.should == 2

        n.shape[1].should == m.shape[0]

        r = n.dot m
        # r.class.should == NVector

        r[0,0].should == 4
        r[1,0].should == 13
        r[2,0].should == 22
        r[3,0].should == 31

        #r.dtype.should == :float64 unless left_dtype == :float32 && right_dtype == :float32
      end
    end
  end


  [:dense, :list, :yale].each do |storage_type|
    context storage_type do
      it "can be duplicated" do
        n = NMatrix.new(storage_type, [2,3], storage_type == :yale ? :float64 : 1.1)
        n.stype.should equal(storage_type)

        n[0,0] = 0.0
        n[0,1] = 0.1
        n[1,0] = 1.0

        m = n.dup
        m.shape.should == n.shape
        m.dim.should == n.dim
        m.object_id.should_not == n.object_id
        m.stype.should equal(storage_type)
        m[0,0].should == n[0,0]
        m[0,0] = 3.0
        m[0,0].should_not == n[0,0]
      end

      it "enforces shape boundaries" do
        lambda { NMatrix.new(storage_type, [1,10], storage_type == :yale ? :int8 : 0)[-1,0] }.should raise_error
        lambda { NMatrix.new(storage_type, [1,10], storage_type == :yale ? :int8 : 0)[1,0]  }.should raise_error(ArgumentError, "out of range")
        lambda { NMatrix.new(storage_type, [1,10], storage_type == :yale ? :int8 : 0)[0,10] }.should raise_error(ArgumentError, "out of range")
      end

      it "sets and gets" do
        n = NMatrix.new(storage_type, 2, storage_type == :yale ? :int8 : 0)
        n[0,1] = 1
        n[0,0].should == 0
        n[1,0].should == 0
        n[0,1].should == 1
        n[1,1].should == 0
      end

      it "sets and gets references" do
        n = NMatrix.new(storage_type, 2, storage_type == :yale ? :int8 : 0)
        (n[0,1] = 1).should == 1
        n[0,1].should == 1
      end

    end

    # dense and list, not yale
    context "(storage: #{storage_type})" do
      it "gets default value" do
        NMatrix.new(storage_type, 3, 0)[1,1].should   == 0
        NMatrix.new(storage_type, 3, 0.1)[1,1].should == 0.1
        NMatrix.new(storage_type, 3, 1)[1,1].should   == 1
      end

      it "returns shape and dim" do
        NMatrix.new(storage_type, [3,2,8], 0).shape.should == [3,2,8]
        NMatrix.new(storage_type, [3,2,8], 0).dim.should  == 3
      end
    end unless storage_type == :yale
  end


  it "handles dense construction" do
    NMatrix.new(3,0)[1,1].should == 0
    lambda { NMatrix.new(3,:int8)[1,1] }.should_not raise_error
  end

  it "allows iteration of Ruby object matrices" do
    n = NMatrix.new(:dense, [3,3], [1,2,3,4,5,6,7,8,9], :object)
    n.each do |x|
      puts x
    end
  end

  it "allows iteration of non-Ruby object matrices" do
    n = NMatrix.new(:dense, [3,3], [1,2,3,4,5,6,7,8,9], :int64)
    n.each do |x|
      puts x
    end
  end

  it "calculates the complex conjugate in-place" do
    n = NMatrix.new(:dense, 3, [1,2,3,4,5,6,7,8,9], :complex128)
    n.complex_conjugate!
    # FIXME: Actually test that values are correct.
  end

  it "converts from list to yale properly" do
    m = NMatrix.new(:list, 3, 0)
    m[0,2] = 333 
    m[2,2] = 777
    n = m.cast(:yale, :int32)
    puts n.capacity
    n.extend NMatrix::YaleFunctions
    puts n.yale_ija.inspect
    puts n.yale_a.inspect
    n[0,0].should == 0
    n[0,1].should == 0
    n[0,2].should == 333
    n[1,0].should == 0
    n[1,1].should == 0
    n[1,2].should == 0
    n[2,0].should == 0
    n[2,1].should == 0
    n[2,2].should == 777
  end

end
