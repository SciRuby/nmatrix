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
# == 00_nmatrix_spec.rb
#
# Basic tests for NMatrix. These should load first, as they're
# essential to NMatrix operation.
#

require File.dirname(__FILE__) + "/spec_helper.rb"

describe NMatrix do
  #after :each do
  #  GC.start
  #end

  it "creates a matrix with the new constructor" do
    n = NMatrix.new([2,2], [0,1,2,3], dtype: :int64)
  end

  it "adequately requires information to access a single entry of a dense matrix" do
    n = NMatrix.new(:dense, 4, [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], :float64)
    n[0,0].should == 0
    expect { n[0] }.to raise_error(ArgumentError)
  end

  it "calculates exact determinants on small square matrices" do
    NMatrix.new(2, [1,2,3,4], stype: :dense, dtype: :int64).det_exact.should == -2
  end

  it "calculates determinants" do
    NMatrix.new(3, [-2,2,3,-1,1,3,2,0,-1], stype: :dense, dtype: :int64).det.should == 6
  end

  it "allows casting to Ruby objects" do
    m = NMatrix.new([3,3], [0,0,1,0,2,0,3,4,5], dtype: :int64, stype: :dense)
    n = m.cast(:dense, :object)
    n.should == m
  end

  it "allows casting from Ruby objects" do
    m = NMatrix.new(:dense, [3,3], [0,0,1,0,2,0,3,4,5], :object)
    n = m.cast(:dense, :int64)
    m.should == n
  end

  it "allows stype casting of a dim 2 matrix between dense, sparse, and list (different dtypes)" do
    m = NMatrix.new(:dense, [3,3], [0,0,1,0,2,0,3,4,5], :int64).
      cast(:yale, :int32).
      cast(:dense, :float64).
      cast(:list, :object).
      cast(:dense, :int16).
      cast(:list, :int32).
      cast(:yale, :int64) #.
    #cast(:list, :int32).
    #cast(:dense, :int16)
    #m.should.equal?(original)
    # For some reason this causes some weird garbage collector problems when we uncomment these. The above lines won't
    # work at all in IRB, but work fine when run in a regular Ruby session.
  end


  it "fills dense Ruby object matrix with nil" do
    n = NMatrix.new([4,3], dtype: :object)
    n[0,0].should == nil
  end

  it "fills dense with individual assignments" do
    n = NMatrix.new([4,3], dtype: :float64)
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
    m = NMatrix.new([4,3], [14.0, 9.0, 3.0, 2.0, 11.0, 15.0, 0.0, 12.0, 17.0, 5.0, 2.0, 3.0], dtype: :float32)
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


  it "dense handles missing initialization value" do
    n = NMatrix.new(3, dtype: :int8)
    n.stype.should == :dense
    n.dtype.should == :int8

    m = NMatrix.new(4, dtype: :float64)
    m.stype.should == :dense
    m.dtype.should == :float64
  end


  [:dense, :list, :yale].each do |storage_type|
    context storage_type do
    it "can be duplicated" do
        n = NMatrix.new([2,3], 1.1, stype: storage_type, dtype: :float64)
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
        expect { NMatrix.new([1,10], 0, dtype: :int8, stype: storage_type, default: 0)[1,0]  }.to raise_error(RangeError)
        expect { NMatrix.new([1,10], 0, dtype: :int8, stype: storage_type, default: 0)[0,10] }.to raise_error(RangeError)
      end

      it "sets and gets" do
        n = NMatrix.new(2, 0, stype: storage_type, dtype: :int8)
        n[0,1] = 1
        n[0,0].should == 0
        n[1,0].should == 0
        n[0,1].should == 1
        n[1,1].should == 0
      end

      it "sets and gets references" do
        n = NMatrix.new(2, stype: storage_type, dtype: :int8, default: 0)
        (n[0,1] = 1).should == 1
        n[0,1].should == 1
      end

      # Tests Ruby object versus any C dtype (in this case we use :int64)
      [:object, :int64].each do |dtype|
        c = dtype == :object ? "Ruby object" : "non-Ruby object"
        context c do
          it "allows iteration of matrices" do
            n = nil
            if storage_type == :dense
              n = NMatrix.new(:dense, [3,3], [1,2,3,4,5,6,7,8,9], dtype)
            else
              n = NMatrix.new([3,4], 0, stype: storage_type, dtype: dtype)
              n[0,0] = 1
              n[0,1] = 2
              n[2,3] = 4
              n[2,0] = 3
            end

            ary = []
            n.each do |x|
              ary << x
            end

            if storage_type == :dense
              ary.should == [1,2,3,4,5,6,7,8,9]
            else
              ary.should == [1,2,0,0,0,0,0,0,3,0,0,4]
            end
          end

          it "allows storage-based iteration of matrices" do
            STDERR.puts storage_type.inspect
            STDERR.puts dtype.inspect
            n = NMatrix.new([3,3], 0, stype: storage_type, dtype: dtype)
            n[0,0] = 1
            n[0,1] = 2
            n[2,0] = 5 if storage_type == :yale
            n[2,1] = 4
            n[2,2] = 3

            values = []
            is = []
            js = []

            n.each_stored_with_indices do |v,i,j|
              values << v
              is << i
              js << j
            end


            if storage_type == :yale
              is.should     == [0,1,2,0,2,2]
              js.should     == [0,1,2,1,0,1]
              values.should == [1,0,3,2,5,4]
            elsif storage_type == :list
              values.should == [1,2,4,3]
              is.should     == [0,0,2,2]
              js.should     == [0,1,1,2]
            elsif storage_type == :dense
              values.should == [1,2,0,0,0,0,0,4,3]
              is.should     == [0,0,0,1,1,1,2,2,2]
              js.should     == [0,1,2,0,1,2,0,1,2]
            end
          end
        end
      end

    end

    # dense and list, not yale
    context "(storage: #{storage_type})" do
      it "gets default value" do
        NMatrix.new(3, 0, stype: storage_type)[1,1].should   == 0
        NMatrix.new(3, 0.1, stype: storage_type)[1,1].should == 0.1
        NMatrix.new(3, 1, stype: storage_type)[1,1].should   == 1
      end

      it "returns shape and dim" do
        NMatrix.new([3,2,8], 0, stype: storage_type).shape.should == [3,2,8]
        NMatrix.new([3,2,8], 0, stype: storage_type).dim.should  == 3
      end

      it "returns number of rows and columns" do
        NMatrix.new([7, 4], 3, stype: storage_type).rows.should == 7
        NMatrix.new([7, 4], 3, stype: storage_type).cols.should == 4
      end
    end unless storage_type == :yale
  end


  it "handles dense construction" do
    NMatrix.new(3,0)[1,1].should == 0
    lambda { NMatrix.new(3,dtype: :int8)[1,1] }.should_not raise_error
  end

  it "calculates the complex conjugate in-place" do
    n = NMatrix.new(:dense, 3, [1,2,3,4,5,6,7,8,9], :complex128)
    n.complex_conjugate!
    # FIXME: Actually test that values are correct.
  end

  it "converts from list to yale properly" do
    m = NMatrix.new(3, 0, stype: :list)
    m[0,2] = 333
    m[2,2] = 777
    n = m.cast(:yale, :int32)
    #puts n.capacity
    #n.extend NMatrix::YaleFunctions
    #puts n.yale_ija.inspect
    #puts n.yale_a.inspect

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

  it "should return an enumerator when each is called without a block" do
    a = NMatrix.new(2, 1)
    b = NMatrix.new(2, [-1,0,1,0])
    enums = [a.each, b.each]

    begin
      atans = []
      atans << Math.atan2(*enums.map(&:next)) while true
    rescue StopIteration
    end
  end

  context "dense" do
    it "should return the matrix being iterated over when each is called with a block" do
      a = NMatrix.new(2, 1)
      val = (a.each { })
      val.should eq a
    end

    it "should return the matrix being iterated over when each_stored_with_indices is called with a block" do
      a = NMatrix.new(2,1)
      val = (a.each_stored_with_indices { })
      val.should eq a
    end
  end

  [:list, :yale].each do |storage_type|
    context storage_type do
      it "should return the matrix being iterated over when each_stored_with_indices is called with a block" do
        n = NMatrix.new([2,3], 1.1, stype: storage_type, dtype: :float64, default: 0)
        val = (n.each_stored_with_indices { })
        val.should eq n
      end

      it "should return an enumerator when each_stored_with_indices is called without a block" do
        n = NMatrix.new([2,3], 1.1, stype: storage_type, dtype: :float64, default: 0)
        val = n.each_stored_with_indices
        val.should be_a Enumerator
      end

    end
  end

  it "should iterate through element 256 without a segfault" do
    t = NVector.random(256)
    t.each { |x| x + 0 }
  end

end


describe "NMatrix#upper_triangle" do
  it "should create a copy with the lower corner set to zero" do
    n = NMatrix.seq(4)+1
    n.upper_triangle.should == NMatrix.new(4, [1,2,3,4,0,6,7,8,0,0,11,12,0,0,0,16])
    n.upper_triangle(2).should == NMatrix.new(4, [1,2,3,4,5,6,7,8,9,10,11,12,0,14,15,16])
  end
end

describe "NMatrix#lower_triangle" do
  it "should create a copy with the lower corner set to zero" do
    n = NMatrix.seq(4)+1
    n.lower_triangle.should == NMatrix.new(4, [1,0,0,0,5,6,0,0,9,10,11,0,13,14,15,16])
    n.lower_triangle(2).should == NMatrix.new(4, [1,2,3,0,5,6,7,8,9,10,11,12,13,14,15,16])
  end
end

describe "NMatrix#upper_triangle!" do
  it "should create a copy with the lower corner set to zero" do
    n = NMatrix.seq(4)+1
    n.upper_triangle!.should == NMatrix.new(4, [1,2,3,4,0,6,7,8,0,0,11,12,0,0,0,16])
    n = NMatrix.seq(4)+1
    n.upper_triangle!(2).should == NMatrix.new(4, [1,2,3,4,5,6,7,8,9,10,11,12,0,14,15,16])
  end
end

describe "NMatrix#lower_triangle!" do
  it "should create a copy with the lower corner set to zero" do
    n = NMatrix.seq(4)+1
    n.lower_triangle!.should == NMatrix.new(4, [1,0,0,0,5,6,0,0,9,10,11,0,13,14,15,16])
    n = NMatrix.seq(4)+1
    n.lower_triangle!(2).should == NMatrix.new(4, [1,2,3,0,5,6,7,8,9,10,11,12,13,14,15,16])
  end
end

describe "NMatrix#reshape" do
  it "should change the shape of a matrix without the contents changing" do
    n = NMatrix.seq(4)+1
    n.reshape([8,2]).to_flat_array.should == n.to_flat_array
  end

  it "should permit a change of dimensionality" do
    n = NMatrix.seq(4)+1
    n.reshape([8,1,2]).to_flat_array.should == n.to_flat_array
  end

  it "should prevent a resize" do
    n = NMatrix.seq(4)+1
    expect { n.reshape([5,2]) }.to raise_error(ArgumentError)
  end
end

describe "NMatrix#transpose" do
  [:dense, :list, :yale].each do |stype|
    context(stype) do
      it "should transpose a #{stype} matrix (2-dimensional)" do
        n = NMatrix.seq(4, stype: stype)
        n.transpose.to_a.flatten.should == [0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15]
      end
    end
  end

  [:dense, :list].each do |stype|
    context(stype) do
      it "should transpose a #{stype} matrix (3-dimensional)" do
        n = NMatrix.new([4,4,1], [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], stype: stype)
        n.transpose([2,1,0]).to_flat_array.should == [0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15]
        n.transpose([1,0,2]).to_flat_array.should == [0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15]
        n.transpose([0,2,1]).to_flat_array.should == n.to_flat_array # for dense, make this reshape!
      end
    end
  end

end

describe "NMatrix#==" do
  [:dense, :list, :yale].each do |left|
    [:dense, :list, :yale].each do |right|
      next if left == right
      context ("#{left}?#{right}") do
        it "should compare two matrices of differing stypes" do
          n = NMatrix.new([3,4], [0,0,1,2,0,0,3,4,0,0,0,0,5,6,7,0], stype: left)
          m = NMatrix.new([3,4], [0,0,1,2,0,0,3,4,0,0,0,0,5,6,7,0], stype: right)
          n.should == m
        end
      end
    end
  end
end

describe "NMatrix#concat" do
  it "should default to horizontal concatenation" do
    n = NMatrix.new([1,3], [1,2,3])
    n.concat(n).should == NMatrix.new([1,6], [1,2,3,1,2,3])
  end

  it "should permit vertical concatenation" do
    n = NMatrix.new([1,3], [1,2,3])
    n.vconcat(n).should == NMatrix.new([2,3], [1,2,3])
  end

  it "should permit depth concatenation on tensors" do
    n = NMatrix.new([1,3,1], [1,2,3])
    n.dconcat(n).should == NMatrix.new([1,3,2], [1,1,2,2,3,3])
  end
end

describe "NMatrix#indexing" do
  it "should return values based on indices" do
    n = NMatrix.new([2,5], [1,2,3,4,5,6,7,8,9,0])
    n[1,0].should == 6
    n[1,0..3].should == NMatrix.new([1,4],[6,7,8,9])
  end

  it "should work for negative indices" do
    n = NMatrix.new([1,5], [1,2,3,4,5])
    m = NMatrix.new([3,3], [1,2,3,4,5,6,7,8,9])
    n[-1].should == 5
    n[0,0..-2].should == NMatrix.new([1,4],[1,2,3,4])
  end
end

