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
require_relative 'spec_helper'

describe NMatrix do
  it "creates a matrix with the new constructor" do
    n = NMatrix.new([2,2], [0,1,2,3], dtype: :int64)
    expect(n.shape).to eq([2,2])
    expect(n.entries).to eq([0,1,2,3])
    expect(n.dtype).to eq(:int64)
  end

  it "adequately requires information to access a single entry of a dense matrix" do
    n = NMatrix.new(:dense, 4, [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], :float64)
    expect(n[0,0]).to eq(0)
    expect { n[0] }.to raise_error(ArgumentError)
  end

  it "calculates exact determinants on small square matrices" do
    expect(NMatrix.new(2, [1,2,3,4], stype: :dense, dtype: :int64).det_exact).to eq(-2)
  end

  it "calculates determinants" do
    expect(NMatrix.new(3, [-2,2,3,-1,1,3,2,0,-1], stype: :dense, dtype: :int64).det).to eq(6)
  end

  it "allows casting to Ruby objects" do
    m = NMatrix.new([3,3], [0,0,1,0,2,0,3,4,5], dtype: :int64, stype: :dense)
    n = m.cast(:dense, :object)
    expect(n).to eq(m)
  end

  it "allows casting from Ruby objects" do
    m = NMatrix.new(:dense, [3,3], [0,0,1,0,2,0,3,4,5], :object)
    n = m.cast(:dense, :int64)
    expect(m).to eq(n)
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
    if jruby?
      expect(n[0,0]).to eq(0.0)
    else
      expect(n[0,0]).to eq(nil)
    end
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

    expect(n[0,0]).to eq(14.0)
    expect(n[0,1]).to eq(9.0)
    expect(n[0,2]).to eq(3.0)
    expect(n[1,0]).to eq(2.0)
    expect(n[1,1]).to eq(11.0)
    expect(n[1,2]).to eq(15.0)
    expect(n[2,0]).to eq(0.0)
    expect(n[2,1]).to eq(12.0)
    expect(n[2,2]).to eq(17.0)
    expect(n[3,0]).to eq(5.0)
    expect(n[3,1]).to eq(2.0)
    expect(n[3,2]).to eq(3.0)
  end

  it "fills dense with a single mass assignment" do
    n = NMatrix.new([4,3], [14.0, 9.0, 3.0, 2.0, 11.0, 15.0, 0.0, 12.0, 17.0, 5.0, 2.0, 3.0])

    expect(n[0,0]).to eq(14.0)
    expect(n[0,1]).to eq(9.0)
    expect(n[0,2]).to eq(3.0)
    expect(n[1,0]).to eq(2.0)
    expect(n[1,1]).to eq(11.0)
    expect(n[1,2]).to eq(15.0)
    expect(n[2,0]).to eq(0.0)
    expect(n[2,1]).to eq(12.0)
    expect(n[2,2]).to eq(17.0)
    expect(n[3,0]).to eq(5.0)
    expect(n[3,1]).to eq(2.0)
    expect(n[3,2]).to eq(3.0)
  end

  it "fills dense with a single mass assignment, with dtype specified" do
    m = NMatrix.new([4,3], [14.0, 9.0, 3.0, 2.0, 11.0, 15.0, 0.0, 12.0, 17.0, 5.0, 2.0, 3.0], dtype: :float32)

    expect(m[0,0]).to eq(14.0)
    expect(m[0,1]).to eq(9.0)
    expect(m[0,2]).to eq(3.0)
    expect(m[1,0]).to eq(2.0)
    expect(m[1,1]).to eq(11.0)
    expect(m[1,2]).to eq(15.0)
    expect(m[2,0]).to eq(0.0)
    expect(m[2,1]).to eq(12.0)
    expect(m[2,2]).to eq(17.0)
    expect(m[3,0]).to eq(5.0)
    expect(m[3,1]).to eq(2.0)
    expect(m[3,2]).to eq(3.0)
  end

  it "dense handles missing initialization value" do
    n = NMatrix.new(3, dtype: :int8)
    expect(n.stype).to eq(:dense)
    expect(n.dtype).to eq(:int8)

    m = NMatrix.new(4, dtype: :float64)
    expect(m.stype).to eq(:dense)
    expect(m.dtype).to eq(:float64)
  end

  [:dense, :list, :yale].each do |storage_type|
    context storage_type do
    it "can be duplicated" do
        n = NMatrix.new([2,3], 1.1, stype: storage_type, dtype: :float64)
        expect(n.stype).to eq(storage_type)

        n[0,0] = 0.0
        n[0,1] = 0.1
        n[1,0] = 1.0

        m = n.dup
        expect(m.shape).to eq(n.shape)
        expect(m.dim).to eq(n.dim)
        expect(m.object_id).not_to eq(n.object_id)
        expect(m.stype).to eq(storage_type)
        expect(m[0,0]).to eq(n[0,0])
        m[0,0] = 3.0
        expect(m[0,0]).not_to eq(n[0,0])
      end

      it "enforces shape boundaries" do
        expect { NMatrix.new([1,10], 0, dtype: :int8, stype: storage_type, default: 0)[1,0]  }.to raise_error(RangeError)
        expect { NMatrix.new([1,10], 0, dtype: :int8, stype: storage_type, default: 0)[0,10] }.to raise_error(RangeError)
      end

      it "sets and gets" do
        n = NMatrix.new(2, 0, stype: storage_type, dtype: :int8)
        n[0,1] = 1
        expect(n[0,0]).to eq(0)
        expect(n[1,0]).to eq(0)
        expect(n[0,1]).to eq(1)
        expect(n[1,1]).to eq(0)
      end

      it "sets and gets references" do
        n = NMatrix.new(2, stype: storage_type, dtype: :int8, default: 0)
        expect(n[0,1] = 1).to eq(1)
        expect(n[0,1]).to eq(1)
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
              expect(ary).to eq([1,2,3,4,5,6,7,8,9])
            else
              expect(ary).to eq([1,2,0,0,0,0,0,0,3,0,0,4])
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
              expect(is).to     eq([0,1,2,0,2,2])
              expect(js).to     eq([0,1,2,1,0,1])
              expect(values).to eq([1,0,3,2,5,4])
            elsif storage_type == :list
              expect(values).to eq([1,2,4,3])
              expect(is).to     eq([0,0,2,2])
              expect(js).to     eq([0,1,1,2])
            elsif storage_type == :dense
              expect(values).to eq([1,2,0,0,0,0,0,4,3])
              expect(is).to     eq([0,0,0,1,1,1,2,2,2])
              expect(js).to     eq([0,1,2,0,1,2,0,1,2])
            end
          end
        end
      end
    end

    # dense and list, not yale
    context "(storage: #{storage_type})" do
      it "gets default value" do
        expect(NMatrix.new(3, 0, stype: storage_type)[1,1]).to eq(0)
        expect(NMatrix.new(3, 0.1, stype: storage_type)[1,1]).to eq(0.1)
        expect(NMatrix.new(3, 1, stype: storage_type)[1,1]).to eq(1)

      end
      it "returns shape and dim" do
        expect(NMatrix.new([3,2,8], 0, stype: storage_type).shape).to eq([3,2,8])
        expect(NMatrix.new([3,2,8], 0, stype: storage_type).dim).to eq(3)
      end

      it "returns number of rows and columns" do
        expect(NMatrix.new([7, 4], 3, stype: storage_type).rows).to eq(7)
        expect(NMatrix.new([7, 4], 3, stype: storage_type).cols).to eq(4)
      end
    end unless storage_type == :yale
  end


  it "handles dense construction" do
    expect(NMatrix.new(3,0)[1,1]).to eq(0)
    expect(lambda { NMatrix.new(3,dtype: :int8)[1,1] }).to_not raise_error
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

    expect(n[0,0]).to eq(0)
    expect(n[0,1]).to eq(0)
    expect(n[0,2]).to eq(333)
    expect(n[1,0]).to eq(0)
    expect(n[1,1]).to eq(0)
    expect(n[1,2]).to eq(0)
    expect(n[2,0]).to eq(0)
    expect(n[2,1]).to eq(0)
    expect(n[2,2]).to eq(777)
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
      expect(val).to eq(a)
    end

    it "should return the matrix being iterated over when each_stored_with_indices is called with a block" do
      a = NMatrix.new(2,1)
      val = (a.each_stored_with_indices { })
      expect(val).to eq(a)
    end
  end

  [:list, :yale].each do |storage_type|
    context storage_type do
      it "should return the matrix being iterated over when each_stored_with_indices is called with a block" do
        n = NMatrix.new([2,3], 1.1, stype: storage_type, dtype: :float64, default: 0)
        val = (n.each_stored_with_indices { })
        expect(val).to eq(n)
      end

      it "should return an enumerator when each_stored_with_indices is called without a block" do
        n = NMatrix.new([2,3], 1.1, stype: storage_type, dtype: :float64, default: 0)
        val = n.each_stored_with_indices
        expect(val).to be_a Enumerator
      end
    end
  end

  it "should iterate through element 256 without a segfault" do
    t = NVector.random(256)
    t.each { |x| x + 0 }
  end
end


describe 'NMatrix' do
  context "#upper_triangle" do
    it "should create a copy with the lower corner set to zero" do
      n = NMatrix.seq(4)+1
      expect(n.upper_triangle).to eq(NMatrix.new(4, [1,2,3,4,0,6,7,8,0,0,11,12,0,0,0,16]))
      expect(n.upper_triangle(2)).to eq(NMatrix.new(4, [1,2,3,4,5,6,7,8,9,10,11,12,0,14,15,16]))
    end
  end

  context "#lower_triangle" do
    it "should create a copy with the lower corner set to zero" do
      n = NMatrix.seq(4)+1
      expect(n.lower_triangle).to eq(NMatrix.new(4, [1,0,0,0,5,6,0,0,9,10,11,0,13,14,15,16]))
      expect(n.lower_triangle(2)).to eq(NMatrix.new(4, [1,2,3,0,5,6,7,8,9,10,11,12,13,14,15,16]))
    end
  end

  context "#upper_triangle!" do
    it "should create a copy with the lower corner set to zero" do
      n = NMatrix.seq(4)+1
      expect(n.upper_triangle!).to eq(NMatrix.new(4, [1,2,3,4,0,6,7,8,0,0,11,12,0,0,0,16]))
      n = NMatrix.seq(4)+1
      expect(n.upper_triangle!(2)).to eq(NMatrix.new(4, [1,2,3,4,5,6,7,8,9,10,11,12,0,14,15,16]))
    end
  end

  context "#lower_triangle!" do
    it "should create a copy with the lower corner set to zero" do
      n = NMatrix.seq(4)+1
      expect(n.lower_triangle!).to eq(NMatrix.new(4, [1,0,0,0,5,6,0,0,9,10,11,0,13,14,15,16]))
      n = NMatrix.seq(4)+1
      expect(n.lower_triangle!(2)).to eq(NMatrix.new(4, [1,2,3,0,5,6,7,8,9,10,11,12,13,14,15,16]))
    end
  end

  context "#rank" do
    it "should get the rank of a 2-dimensional matrix" do
      n = NMatrix.seq([2,3])
      expect(n.rank(0, 0)).to eq(N[[0,1,2]])
    end

    it "should raise an error when the rank is out of bounds" do
      n = NMatrix.seq([2,3])
      expect { n.rank(2, 0) }.to raise_error(RangeError)
    end
  end

  context "#reshape" do
    it "should change the shape of a matrix without the contents changing" do
      n = NMatrix.seq(4)+1
      expect(n.reshape([8,2]).to_flat_array).to eq(n.to_flat_array)
    end

    it "should permit a change of dimensionality" do
      n = NMatrix.seq(4)+1
      expect(n.reshape([8,1,2]).to_flat_array).to eq(n.to_flat_array)
    end

    it "should prevent a resize" do
      n = NMatrix.seq(4)+1
      expect { n.reshape([5,2]) }.to raise_error(ArgumentError)
    end

    it "should do the reshape operation in place" do
      n = NMatrix.seq(4)+1
      expect(n.reshape!([8,2]).eql?(n)).to eq(true) # because n itself changes
    end

    it "should do the reshape operation in place, changing dimension" do
      n = NMatrix.seq(4)
      a = n.reshape!([4,2,2])
      expect(n).to eq(NMatrix.seq([4,2,2]))
      expect(a).to eq(NMatrix.seq([4,2,2]))
    end

    it "reshape and reshape! must produce same result" do
      n = NMatrix.seq(4)+1
      a = NMatrix.seq(4)+1
      expect(n.reshape!([8,2])==a.reshape(8,2)).to eq(true) # because n itself changes
    end

    it "should prevent a resize in place" do
      n = NMatrix.seq(4)+1
      expect { n.reshape!([5,2]) }.to raise_error(ArgumentError)
    end
  end

  context "#transpose" do
    [:dense, :list, :yale].each do |stype|
      context(stype) do
        it "should transpose a #{stype} matrix (2-dimensional)" do
          n = NMatrix.seq(4, stype: stype)
          expect(n.transpose.to_a.flatten).to eq([0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15])
        end
      end
    end

    [:dense, :list].each do |stype|
      context(stype) do
        it "should transpose a #{stype} matrix (3-dimensional)" do
          n = NMatrix.new([4,4,1], [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], stype: stype)
          expect(n.transpose([2,1,0]).to_flat_array).to eq([0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15])
          expect(n.transpose([1,0,2]).to_flat_array).to eq([0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15])
          expect(n.transpose([0,2,1]).to_flat_array).to eq(n.to_flat_array) # for dense, make this reshape!
        end
      end

      it "should just copy a 1-dimensional #{stype} matrix" do
        n = NMatrix.new([3], [1,2,3], stype: stype)
        expect(n.transpose).to eq n
        expect(n.transpose).not_to be n
      end

      it "should check permute argument if supplied for #{stype} matrix" do
        n = NMatrix.new([2,2], [1,2,3,4], stype: stype)
        expect{n.transpose *4 }.to raise_error(ArgumentError)
        expect{n.transpose [1,1,2] }.to raise_error(ArgumentError)
      end
    end
  end

  context "#dot_product" do
    [:dense].each do |stype| # list storage transpose not yet implemented
      context(stype) do # yale support only 2-dim matrix
        it "should work like vector product on a #{stype} (1-dimensional)" do
          m = NMatrix.new([3], [1,2,3], stype: stype)
          expect(m.dot(m)).to eq (NMatrix.new([1],[14]))
        end
      end
    end
  end

  context "#==" do
    [:dense, :list, :yale].each do |left|
      [:dense, :list, :yale].each do |right|
        context ("#{left}?#{right}") do
          it "tests equality of two equal matrices" do
            n = NMatrix.new([3,4], [0,0,1,2,0,0,3,4,0,0,0,0], stype: left)
            m = NMatrix.new([3,4], [0,0,1,2,0,0,3,4,0,0,0,0], stype: right)

            expect(n==m).to eq(true)
          end

          it "tests equality of two unequal matrices" do
            n = NMatrix.new([3,4], [0,0,1,2,0,0,3,4,0,0,0,1], stype: left)
            m = NMatrix.new([3,4], [0,0,1,2,0,0,3,4,0,0,0,0], stype: right)

            expect(n==m).to eq(false)
          end

          it "tests equality of matrices with different shapes" do
            n = NMatrix.new([2,2], [1,2, 3,4], stype: left)
            m = NMatrix.new([2,3], [1,2, 3,4, 5,6], stype: right)
            x = NMatrix.new([1,4], [1,2, 3,4], stype: right)

            expect{n==m}.to raise_error(ShapeError)
            expect{n==x}.to raise_error(ShapeError)
          end

          it "tests equality of matrices with different dimension" do
            n = NMatrix.new([2,1], [1,2], stype: left)
            m = NMatrix.new([2], [1,2], stype: right)

            expect{n==m}.to raise_error(ShapeError)
          end if left != :yale && right != :yale # yale must have dimension 2
        end
      end
    end
  end

  context "#concat" do
    it "should default to horizontal concatenation" do
      n = NMatrix.new([1,3], [1,2,3])
      expect(n.concat(n)).to eq(NMatrix.new([1,6], [1,2,3,1,2,3]))
    end

    it "should permit vertical concatenation" do
      n = NMatrix.new([1,3], [1,2,3])
      expect(n.vconcat(n)).to eq(NMatrix.new([2,3], [1,2,3]))
    end

    it "should permit depth concatenation on tensors" do
      n = NMatrix.new([1,3,1], [1,2,3])
      expect(n.dconcat(n)).to eq(NMatrix.new([1,3,2], [1,1,2,2,3,3]))
    end

    it "should work on matrices with different size along concat dim" do
      n = N[[1, 2, 3],
            [4, 5, 6]]
      m = N[[7],
            [8]]

      expect(n.hconcat(m)).to eq N[[1, 2, 3, 7], [4, 5, 6, 8]]
      expect(m.hconcat(n)).to eq N[[7, 1, 2, 3], [8, 4, 5, 6]]
    end

    it "should work on matrices with different size along concat dim" do
      n = N[[1, 2, 3],
            [4, 5, 6]]

      m = N[[7, 8, 9]]

      expect(n.vconcat(m)).to eq N[[1, 2, 3], [4, 5, 6], [7, 8, 9]]
      expect(m.vconcat(n)).to eq N[[7, 8, 9], [1, 2, 3], [4, 5, 6]]
    end
  end

  context "#[]" do
    it "should return values based on indices" do
      n = NMatrix.new([2,5], [1,2,3,4,5,6,7,8,9,0])
      expect(n[1,0]).to eq 6
      expect(n[1,0..3]).to eq NMatrix.new([1,4],[6,7,8,9])
    end

    it "should work for negative indices" do
      n = NMatrix.new([1,5], [1,2,3,4,5])
      expect(n[-1]).to eq(5)
      expect(n[0,0..-2]).to eq(NMatrix.new([1,4],[1,2,3,4]))
    end
  end

  context "#complex_conjugate!" do
    [:dense, :yale, :list].each do |stype|
      context(stype) do
        it "should work in-place for complex dtypes" do
          pending("not yet implemented for list stype") if stype == :list
          n = NMatrix.new([2,3], [Complex(2,3)], stype: stype, dtype: :complex128)
          n.complex_conjugate!
          expect(n).to eq(NMatrix.new([2,3], [Complex(2,-3)], stype: stype, dtype: :complex128))
        end

        [:object, :int64].each do |dtype|
          it "should work in-place for non-complex dtypes" do
            pending("not yet implemented for list stype") if stype == :list
            n = NMatrix.new([2,3], 1, stype: stype, dtype: dtype)
            n.complex_conjugate!
            expect(n).to eq(NMatrix.new([2,3], [1], stype: stype, dtype: dtype))
          end
        end
      end
    end
  end

  context "#complex_conjugate" do
    [:dense, :yale, :list].each do |stype|
      context(stype) do
        it "should work out-of-place for complex dtypes" do
          pending("not yet implemented for list stype") if stype == :list
          n = NMatrix.new([2,3], [Complex(2,3)], stype: stype, dtype: :complex128)
          expect(n.complex_conjugate).to eq(NMatrix.new([2,3], [Complex(2,-3)], stype: stype, dtype: :complex128))
        end

        [:object, :int64].each do |dtype|
          it "should work out-of-place for non-complex dtypes" do
            pending("not yet implemented for list stype") if stype == :list
            n = NMatrix.new([2,3], 1, stype: stype, dtype: dtype)
            expect(n.complex_conjugate).to eq(NMatrix.new([2,3], [1], stype: stype, dtype: dtype))
          end
        end
      end
    end
  end

  context "#inject" do
    it "should sum columns of yale matrix correctly" do
      n = NMatrix.new([4, 3], stype: :yale, default: 0)
      n[0,0] = 1
      n[1,1] = 2
      n[2,2] = 4
      n[3,2] = 8
      column_sums = []
      n.cols.times do |i|
        column_sums << n.col(i).inject(:+)
      end
      expect(column_sums).to eq([1, 2, 12])
    end
  end

  context "#index" do
    it "returns index of first occurence of an element for a vector" do
      n = NMatrix.new([5], [0,22,22,11,11])

      expect(n.index(22)).to eq([1])
    end

    it "returns index of first occurence of an element for 2-D matrix" do
      n = NMatrix.new([3,3], [23,11,23,
                              44, 2, 0,
                              33, 0, 32])

      expect(n.index(0)).to eq([1,2])
    end

    it "returns index of first occerence of an element for N-D matrix" do
      n = NMatrix.new([3,3,3], [23,11,23, 44, 2, 0, 33, 0, 32,
                                23,11,23, 44, 2, 0, 33, 0, 32,
                                23,11,23, 44, 2, 0, 33, 0, 32])

      expect(n.index(44)).to eq([0,1,0])
    end
  end

  context "#last" do
    it "returns the last element of a 1-dimensional NMatrix" do
      n = NMatrix.new([1,4], [1,2,3,4])
      expect(n.last).to eq(4)
    end

    it "returns the last element of a 2-dimensional NMatrix" do
      n = NMatrix.new([2,2], [4,8,12,16])
      expect(n.last).to eq(16)
    end

    it "returns the last element of a 3-dimensional NMatrix" do
      n = NMatrix.new([2,2,2], [1,2,3,4,5,6,7,8])
      expect(n.last).to eq(8)
    end
  end

  context "#diagonal" do
    ALL_DTYPES.each do |dtype|
      before do 
        @square_matrix =  NMatrix.new([3,3], [
          23,11,23,
          44, 2, 0,
          33, 0, 32
          ], dtype: dtype
        )

        @rect_matrix = NMatrix.new([4,3], [
          23,11,23,
          44, 2, 0,
          33, 0,32,
          11,22,33
          ], dtype: dtype
        )
      end

      it "returns main diagonal for square matrix" do
        expect(@square_matrix.diagonal).to eq(NMatrix.new [3], [23,2,32])
      end

      it "returns main diagonal for rectangular matrix" do
        expect(@rect_matrix.diagonal).to eq(NMatrix.new [3], [23,2,32])
      end

      it "returns anti-diagonal for square matrix" do
        expect(@square_matrix.diagonal(false)).to eq(NMatrix.new [3], [23,2,33])
      end

      it "returns anti-diagonal for rectangular matrix" do
        expect(@square_matrix.diagonal(false)).to eq(NMatrix.new [3], [23,2,33])
      end
    end
  end

  context "#repeat" do
    before do
      @sample_matrix = NMatrix.new([2, 2], [1, 2, 3, 4])
    end

    it "checks count argument" do
      expect{@sample_matrix.repeat(1, 0)}.to raise_error(ArgumentError)
      expect{@sample_matrix.repeat(-2, 0)}.to raise_error(ArgumentError)
    end

    it "returns repeated matrix" do
      expect(@sample_matrix.repeat(2, 0)).to eq(NMatrix.new([4, 2], [1, 2, 3, 4, 1, 2, 3, 4]))
      expect(@sample_matrix.repeat(2, 1)).to eq(NMatrix.new([2, 4], [1, 2, 1, 2, 3, 4, 3, 4]))
    end

    it "preserves dtype" do
      expect(@sample_matrix.repeat(2, 0).dtype).to eq(@sample_matrix.dtype)
      expect(@sample_matrix.repeat(2, 1).dtype).to eq(@sample_matrix.dtype)
    end
  end

  context "#meshgrid" do
    before do
      @x, @y, @z = [1, 2, 3], NMatrix.new([2, 1], [4, 5]), [6, 7]
      @two_dim = NMatrix.new([2, 2], [1, 2, 3, 4])
      @two_dim_array = [[4], [5]]
      @expected_result = [NMatrix.new([2, 3], [1, 2, 3, 1, 2, 3]), NMatrix.new([2, 3], [4, 4, 4, 5, 5, 5])]
      @expected_for_ij = [NMatrix.new([3, 2], [1, 1, 2, 2, 3, 3]), NMatrix.new([3, 2], [4, 5, 4, 5, 4, 5])]
      @expected_for_sparse = [NMatrix.new([1, 3], [1, 2, 3]), NMatrix.new([2, 1], [4, 5])]
      @expected_for_sparse_ij = [NMatrix.new([3, 1], [1, 2, 3]), NMatrix.new([1, 2], [4, 5])]
      @expected_3dim = [NMatrix.new([1, 3, 1], [1, 2, 3]).repeat(2, 0).repeat(2, 2),
                        NMatrix.new([2, 1, 1], [4, 5]).repeat(3, 1).repeat(2, 2),
                        NMatrix.new([1, 1, 2], [6, 7]).repeat(2, 0).repeat(3, 1)]
      @expected_3dim_sparse_ij = [NMatrix.new([3, 1, 1], [1, 2, 3]),
                                  NMatrix.new([1, 2, 1], [4, 5]),
                                  NMatrix.new([1, 1, 2], [6, 7])]
    end

    it "checks arrays count" do
      expect{NMatrix.meshgrid([@x])}.to raise_error(ArgumentError)
      expect{NMatrix.meshgrid([])}.to raise_error(ArgumentError)
    end

    it "flattens input arrays before use" do
      expect(NMatrix.meshgrid([@two_dim, @two_dim_array])).to eq(NMatrix.meshgrid([@two_dim.to_flat_array, @two_dim_array.flatten]))
    end

    it "returns new NMatrixes" do
      expect(NMatrix.meshgrid([@x, @y])).to eq(@expected_result)
    end

    it "has option :sparse" do
      expect(NMatrix.meshgrid([@x, @y], sparse: true)).to eq(@expected_for_sparse)
    end

    it "has option :indexing" do
      expect(NMatrix.meshgrid([@x, @y], indexing: :ij)).to eq(@expected_for_ij)
      expect(NMatrix.meshgrid([@x, @y], indexing: :xy)).to eq(@expected_result)
      expect{NMatrix.meshgrid([@x, @y], indexing: :not_ij_not_xy)}.to raise_error(ArgumentError)
    end

    it "works well with both options set" do
      expect(NMatrix.meshgrid([@x, @y], sparse: true, indexing: :ij)).to eq(@expected_for_sparse_ij)
    end

    it "is able to take more than two arrays as arguments and works well with options" do
      expect(NMatrix.meshgrid([@x, @y, @z])).to eq(@expected_3dim)
      expect(NMatrix.meshgrid([@x, @y, @z], sparse: true, indexing: :ij)).to eq(@expected_3dim_sparse_ij)
    end
  end
end
