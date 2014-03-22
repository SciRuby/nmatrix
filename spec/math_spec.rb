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
# == math_spec.rb
#
# Tests for non-BLAS and non-LAPACK math functions, or for simplified
# versions of unfriendly BLAS and LAPACK functions.
#

require 'spec_helper'

describe "math" do
  #after :each do
  #  GC.start
  #end

  context "elementwise math functions" do

    [:dense,:list,:yale].each do |stype|
      context stype do

        [:int64,:float64,:rational128].each do |dtype|
          context dtype do
            before :each do
              @size = [2,2]
              @m = NMatrix.seq(@size, dtype: dtype, stype: stype)+1
              @a = @m.to_a.flatten
            end

            NMatrix::NMMath::METHODS_ARITY_1.each do |meth|
              #skip inverse regular trig functions
              next if meth.to_s.start_with?('a') and (not meth.to_s.end_with?('h')) \
                and NMatrix::NMMath::METHODS_ARITY_1.include?(
                  meth.to_s[1...meth.to_s.length].to_sym)
              next if meth == :atanh

              it "should correctly apply elementwise #{meth}" do

                expect(@m.send(meth)).to eq N.new(@size, @a.map{ |e| Math.send(meth, e) },
                                                 dtype: :float64, stype: stype)
              end
            end

            NMatrix::NMMath::METHODS_ARITY_2.each do |meth|
              next if meth == :atan2
              it "should correctly apply elementwise #{meth}" do
                expect(@m.send(meth, @m)).to eq N.new(@size, @a.map{ |e|
                                                     Math.send(meth, e, e) },
                                                     dtype: :float64, 
                                                     stype: stype)
              end

              it "should correctly apply elementwise #{meth} with a scalar first arg" do
                expect(Math.send(meth, 1, @m)).to eq N.new(@size, @a.map { |e| Math.send(meth, 1, e) }, dtype: :float64, stype: stype)
              end

              it "should correctly apply elementwise #{meth} with a scalar second arg" do
                expect(@m.send(meth, 1)).to eq N.new(@size, @a.map { |e| Math.send(meth, e, 1) }, dtype: :float64, stype: stype)
              end
            end

            it "should correctly apply elementwise natural log" do
              expect(@m.log).to eq N.new(@size, [0, Math.log(2), Math.log(3), Math.log(4)],
                                        dtype: :float64, stype: stype)
            end

            it "should correctly apply elementwise log with arbitrary base" do
              expect(@m.log(3)).to eq N.new(@size, [0, Math.log(2,3), 1, Math.log(4,3)],
                                           dtype: :float64, stype: stype)
            end

            context "inverse trig functions" do
              before :each do
                @m = NMatrix.seq(@size, dtype: dtype, stype: stype)/4
                @a = @m.to_a.flatten
              end
              [:asin, :acos, :atan, :atanh].each do |atf|

                it "should correctly apply elementwise #{atf}" do
                  expect(@m.send(atf)).to eq N.new(@size, 
                                               @a.map{ |e| Math.send(atf, e) },
                                               dtype: :float64, stype: stype)
                end
              end

              it "should correctly apply elementtwise atan2" do
                expect(@m.atan2(@m*0+1)).to eq N.new(@size, 
                  @a.map { |e| Math.send(:atan2, e, 1) }, dtype: :float64, stype: stype)
              end

              it "should correctly apply elementwise atan2 with a scalar first arg" do
                expect(Math.atan2(1, @m)).to eq N.new(@size, @a.map { |e| Math.send(:atan2, 1, e) }, dtype: :float64, stype: stype)
              end

              it "should correctly apply elementwise atan2 with a scalar second arg" do
                  expect(@m.atan2(1)).to eq N.new(@size, @a.map { |e| Math.send(:atan2, e, 1) }, dtype: :float64, stype: stype)
              end
            end
          end
        end

      end
    end
  end

  [:float32, :float64, :complex64, :complex128, :rational32, :rational64, :rational128].each do |dtype|
    context dtype do
      it "should correctly factorize a matrix" do
        m = NMatrix.new(:dense, 3, [4,9,2,3,5,7,8,1,6], dtype)
        a = m.factorize_lu
        expect(a[0,0]).to eq(8)
        expect(a[0,1]).to eq(1)
        expect(a[0,2]).to eq(6)
        expect(a[1,0]).to eq(0.5)
        expect(a[1,1]).to eq(8.5)
        expect(a[1,2]).to eq(-1)
        expect(a[2,0]).to eq(0.375)
      end
    end

    context dtype do
      it "should correctly invert a matrix in place" do
        a = NMatrix.new(:dense, 3, [1,0,4,1,1,6,-3,0,-10], dtype)
        b = NMatrix.new(:dense, 3, [-5,0,-2,-4,1,-1,3.quo(2),0,1.quo(2)], dtype)
        begin
          a.invert!
        rescue NotImplementedError => e
          if dtype.to_s =~ /rational/
            pending "getri needs rational implementation"
          else
            pending e.to_s
          end
        end
        expect(a).to eq(b)
      end

      unless NMatrix.has_clapack?
        it "should correctly exact-invert a matrix" do
          a = NMatrix.new(:dense, 3, [1,0,4,1,1,6,-3,0,-10], dtype)
          b = NMatrix.new(:dense, 3, [-5,0,-2,-4,1,-1,3.quo(2),0,1.quo(2)], dtype)
          a.invert.should == b
        end
      end
    end
  end
  
  #TODO: check validity of SVD values in 2-norm
  context "#norm" do
    it "should default to 2-norm" do
      n = NMatrix.new([3, 3], [-4, -3, -2, -1, 0, 1, 2, 3, 4])
      expect(n.norm).to eq(7.348469228349535)
    end

    it "should reject invalid arguments" do
      n = NMatrix.new([3, 3], [-4, -3, -2, -1, 0, 1, 2, 3, 4])
                 
      expect{n.norm(0.5)}.to raise_error(ArgumentError)
    end
    
    it "should calculate 1 and 2 norms correctly" do
      n = NMatrix.new([3, 3], [-4, -3, -2, -1, 0, 1, 2, 3, 4])
      expect(n.norm(1)).to eq(7)
      #FIXME: change to the correct value when overflow issue is resolved
      #expect(n.norm(-2)).to eq(1.8628605857884395e-07)
      expect(n.norm(-2)).to eq(0.0)
      expect(n.norm(-1)).to eq(6)
    end
    
    it "should calculate infinity norms correctly" do
      n = NMatrix.new([3, 3], [-4, -3, -2, -1, 0, 1, 2, 3, 4])
      expect(n.norm(:inf)).to eq(9)
      expect(n.norm(:'-inf')).to eq(2)
    end
    
    it "should calculate frobenius norms correctly" do
      n = NMatrix.new([3, 3], [-4, -3, -2, -1, 0, 1, 2, 3, 4])
      expect(n.norm(:fro)).to eq(7.745966692414834)
    end
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

        n = NMatrix.new([4,3], nary, dtype: left_dtype, stype: :dense)
        m = NMatrix.new([3,2], mary, dtype: right_dtype, stype: :dense)

        expect(m.shape[0]).to eq(3)
        expect(m.shape[1]).to eq(2)
        expect(m.dim).to eq(2)

        expect(n.shape[0]).to eq(4)
        expect(n.shape[1]).to eq(3)
        expect(n.dim).to eq(2)

        expect(n.shape[1]).to eq(m.shape[0])

        r = n.dot m

        expect(r[0,0]).to eq(273.0)
        expect(r[0,1]).to eq(455.0)
        expect(r[1,0]).to eq(243.0)
        expect(r[1,1]).to eq(235.0)
        expect(r[2,0]).to eq(244.0)
        expect(r[2,1]).to eq(205.0)
        expect(r[3,0]).to eq(102.0)
        expect(r[3,1]).to eq(160.0)

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
        n = NMatrix.new([4,3], [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0], dtype: left_dtype)

        m = NMatrix.new([3,1], [2.0, 1.0, 0.0], dtype: right_dtype)

        expect(m.shape[0]).to eq(3)
        expect(m.shape[1]).to eq(1)

        expect(n.shape[0]).to eq(4)
        expect(n.shape[1]).to eq(3)
        expect(n.dim).to eq(2)

        expect(n.shape[1]).to eq(m.shape[0])

        r = n.dot m
        # r.class.should == NVector

        expect(r[0,0]).to eq(4)
        expect(r[1,0]).to eq(13)
        expect(r[2,0]).to eq(22)
        expect(r[3,0]).to eq(31)

        #r.dtype.should == :float64 unless left_dtype == :float32 && right_dtype == :float32
      end
    end
  end
end
