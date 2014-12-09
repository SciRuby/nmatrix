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

ALL_DTYPES = [:byte,:int8,:int16,:int32,:int64, :float32,:float64, :object,
  :rational32,:rational64,:rational128, :complex64, :complex128]

NON_INTEGER_DTYPES = [:float32, :float64, :complex64, :complex128, :rational32, :rational64, :rational128, :object]

describe "math" do
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

              if meth == :-@
                it "should correctly apply elementwise negation" do
                  expect(@m.send(meth)).to eq N.new(@size, @a.map { |e| -e }, dtype: dtype, stype: stype)
                end
                next
              end

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
          
        context "Floor and ceil for #{stype}" do  

          [:floor, :ceil].each do |meth|
            ALL_DTYPES.each do |dtype|
              context dtype do
                before :each do
                  @size = [2,2]
                  @m    = NMatrix.seq(@size, dtype: dtype, stype: stype)+1
                  @a    = @m.to_a.flatten
                end

                if dtype.to_s.match(/int/) or [:byte, :object].include?(dtype)
                  it "should return #{dtype} for #{dtype}" do
                    
                    expect(@m.send(meth)).to eq N.new(@size, @a.map { |e| e.send(meth) }, dtype: dtype, stype: stype)

                    if dtype == :object
                      expect(@m.send(meth).dtype).to eq :object
                    else
                      expect(@m.send(meth).integer_dtype?).to eq true
                    end
                  end
                elsif dtype.to_s.match(/float/) or dtype.to_s.match(/rational/) 
                  it "should return dtype int64 for #{dtype}" do

                    expect(@m.send(meth)).to eq N.new(@size, @a.map { |e| e.send(meth) }, dtype: dtype, stype: stype)
                    
                    expect(@m.send(meth).dtype).to eq :int64
                  end
                elsif dtype.to_s.match(/complex/) 
                  it "should properly calculate #{meth} for #{dtype}" do

                    expect(@m.send(meth)).to eq N.new(@size, @a.map { |e| e = Complex(e.real.send(meth), e.imag.send(meth)) }, dtype: dtype, stype: stype)

                    expect(@m.send(meth).dtype).to eq :complex64  if dtype == :complex64
                    expect(@m.send(meth).dtype).to eq :complex128 if dtype == :complex128
                  end
                end
              end
            end
          end
        end

        context "#round for #{stype}" do
          ALL_DTYPES.each do |dtype|
            context dtype do
              before :each do
                @size = [2,2]
                @mat  = NMatrix.new @size, [1.33334, 0.9998, 1.9999, -8.9999], 
                  dtype: dtype, stype: stype
                @ans  = @mat.to_a.flatten
              end

              it "rounds #{dtype} for #{stype}" do
                expect(@mat.round).to eq(N.new(@size, @ans.map { |a| a.round}, 
                  dtype: dtype, stype: stype))
              end unless(/complex/ =~ dtype)

              it "rounds complex dtype #{dtype} for #{stype}" do
                
                expect(@mat.round).to eq(N.new [2,2], @ans.map {|a| 
                  Complex(a.real.round, a.imag.round)},dtype: dtype, stype: stype)
              end if(/complex/ =~ dtype)
            end
          end
        end
        
      end
    end
  end

  NON_INTEGER_DTYPES.each do |dtype|
    next if dtype == :object
    context dtype do
      before do
        @m = NMatrix.new(:dense, 3, [4,9,2,3,5,7,8,1,6], dtype)
      end

      it "should correctly factorize a matrix" do
        a = @m.factorize_lu
        expect(a[0,0]).to eq(8)
        expect(a[0,1]).to eq(1)
        expect(a[0,2]).to eq(6)
        expect(a[1,0]).to eq(0.5)
        expect(a[1,1]).to eq(8.5)
        expect(a[1,2]).to eq(-1)
        expect(a[2,0]).to eq(0.375)
      end

      it "also returns the permutation matrix" do
        a, p = @m.factorize_lu perm_matrix: true

        expect(a[0,0]).to eq(8)
        expect(a[0,1]).to eq(1)
        expect(a[0,2]).to eq(6)
        expect(a[1,0]).to eq(0.5)
        expect(a[1,1]).to eq(8.5)
        expect(a[1,2]).to eq(-1)
        expect(a[2,0]).to eq(0.375)

        puts p
        expect(p[1,0]).to eq(1)
        expect(p[2,1]).to eq(1)
        expect(p[0,2]).to eq(1)
      end
    end

    context dtype do
      it "should correctly invert a matrix in place (bang)", :focus => true do
        a = NMatrix.new(:dense, 3, [1,2,3,0,1,4,5,6,0], dtype)
        b = NMatrix.new(:dense, 3, [-24,18,5,20,-15,-4,-5,4,1], dtype)
        begin
          a.invert!
        rescue NotImplementedError => e
          if dtype.to_s =~ /rational/
            pending "getri needs rational implementation"
          else
            pending e.to_s
          end
        end
        expect(a.round).to eq(b)
      end

      unless NMatrix.has_clapack?
        it "should correctly invert a matrix in place" do
          a = NMatrix.new(:dense, 5, [1, 8,-9, 7, 5, 
                                      0, 1, 0, 4, 4, 
                                      0, 0, 1, 2, 5, 
                                      0, 0, 0, 1,-5,
                                      0, 0, 0, 0, 1 ], dtype)
          b = NMatrix.new(:dense, 5, [1,-8, 9, 7, 17,
                                      0, 1, 0,-4,-24,
                                      0, 0, 1,-2,-15,
                                      0, 0, 0, 1,  5,
                                      0, 0, 0, 0,  1,], dtype)
          expect(a.invert).to eq(b)
        end
      end

      it "should correctly invert a matrix out-of-place" do
        a = NMatrix.new(:dense, 3, [1,2,3,0,1,4,5,6,0], dtype)
        b = NMatrix.new(:dense, 3, [-24,18,5,20,-15,-4,-5,4,1], dtype)

        expect(a.invert(3,3)).to eq(b)
      end
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

  context "#solve" do
    NON_INTEGER_DTYPES.each do |dtype|
      next if dtype == :object # LU factorization doesnt work for :object yet
      it "solves linear equation for dtype #{dtype}" do
        a = NMatrix.new [2,2], [3,1,1,2], dtype: dtype
        b = NMatrix.new [2,1], [9,8], dtype: dtype

        expect(a.solve(b)).to eq(NMatrix.new [2,1], [2,3], dtype: dtype)
      end
    end
  end
end
