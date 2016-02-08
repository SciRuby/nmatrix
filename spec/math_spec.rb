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
  context "elementwise math functions" do

    [:dense,:list,:yale].each do |stype|
      context stype do

        [:int64,:float64].each do |dtype|
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
                elsif dtype.to_s.match(/float/)
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

              it "rounds" do
                expect(@mat.round).to eq(N.new(@size, @ans.map { |a| a.round}, 
                  dtype: dtype, stype: stype))
              end unless(/complex/ =~ dtype)

              it "rounds with args" do
                expect(@mat.round(2)).to eq(N.new(@size, @ans.map { |a| a.round(2)}, 
                  dtype: dtype, stype: stype))
              end unless(/complex/ =~ dtype)

              it "rounds complex with args" do
                puts @mat.round(2)
                expect(@mat.round(2)).to be_within(0.0001).of(N.new [2,2], @ans.map {|a| 
                  Complex(a.real.round(2), a.imag.round(2))},dtype: dtype, stype: stype)
              end if(/complex/ =~ dtype)

              it "rounds complex" do
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
        @m = NMatrix.new([3,4], GETRF_EXAMPLE_ARRAY, dtype: dtype)
        @err = case dtype
                 when :float32, :complex64
                   1e-6
                 when :float64, :complex128
                   1e-14
               end
      end

      #haven't check this spec yet. Also it doesn't check all the elements of the matrix.
      it "should correctly factorize a matrix" do
        a = @m.factorize_lu
        expect(a).to be_within(@err).of(NMatrix.new([3,4], GETRF_SOLUTION_ARRAY, dtype: dtype))
      end

      it "also returns the permutation matrix" do
        a, p = @m.factorize_lu perm_matrix: true

        expect(a).to be_within(@err).of(NMatrix.new([3,4], GETRF_SOLUTION_ARRAY, dtype: dtype))

        p_true = NMatrix.new([3,3], [0,0,1,1,0,0,0,1,0], dtype: dtype)
        expect(p).to eq(p_true)
      end
    end
  end

  NON_INTEGER_DTYPES.each do |dtype|
    next if dtype == :object
    context dtype do

      it "calculates cholesky decomposition using potrf (lower)" do
        #a = NMatrix.new([3,3],[1,1,1, 1,2,2, 1,2,6], dtype: dtype)
        # We use the matrix
        # 1 1 1
        # 1 2 2
        # 1 2 6
        # which is symmetric and positive-definite as required, but
        # we need only store the lower-half of the matrix.
        a = NMatrix.new([3,3],[1,0,0, 1,2,0, 1,2,6], dtype: dtype)
        begin
          r = a.potrf!(:lower)

          b = NMatrix.new([3,3],[1,0,0, 1,1,0, 1,1,2], dtype: dtype)
          expect(a).to eq(b)
          expect(r).to eq(b)
        rescue NotImplementedError
          pending "potrf! not implemented without plugins"
        end
      end

      it "calculates cholesky decomposition using potrf (upper)" do
        a = NMatrix.new([3,3],[1,1,1, 0,2,2, 0,0,6], dtype: dtype)
        begin
          r = a.potrf!(:upper)

          b = NMatrix.new([3,3],[1,1,1, 0,1,1, 0,0,2], dtype: dtype)
          expect(a).to eq(b)
          expect(r).to eq(b)
        rescue NotImplementedError
          pending "potrf! not implemented without plugins"
        end
      end

      it "calculates cholesky decomposition using #factorize_cholesky" do
        a = NMatrix.new([3,3],[1,2,1, 2,13,5, 1,5,6], dtype: dtype)
        begin
          u,l = a.factorize_cholesky

          l_true = NMatrix.new([3,3],[1,0,0, 2,3,0, 1,1,2], dtype: dtype)
          u_true = l_true.transpose
          expect(u).to eq(u_true)
          expect(l).to eq(l_true)
        rescue NotImplementedError
          pending "potrf! not implemented without plugins"
        end
      end
    end
  end

  ALL_DTYPES.each do |dtype|
    next if dtype == :byte #doesn't work for unsigned types
    next if dtype == :object

    context dtype do
      err = case dtype
              when :float32, :complex64
                1e-4
              else #integer matrices will return :float64
                1e-13
            end

      it "should correctly invert a matrix in place (bang)" do
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
        if a.integer_dtype?
          expect{a.invert!}.to raise_error(DataTypeError)
        else
          #should return inverse as well as modifying a
          r = a.invert!
          expect(a).to be_within(err).of(b)
          expect(r).to be_within(err).of(b)
        end
      end

      it "should correctly invert a matrix out-of-place" do
        a = NMatrix.new(:dense, 3, [1,2,3,0,1,4,5,6,0], dtype)

        if a.integer_dtype?
          b = NMatrix.new(:dense, 3, [-24,18,5,20,-15,-4,-5,4,1], :float64)
        else
          b = NMatrix.new(:dense, 3, [-24,18,5,20,-15,-4,-5,4,1], dtype)
        end

        expect(a.invert).to be_within(err).of(b)
      end
    end
  end

  # TODO: Get it working with ROBJ too
  [:byte,:int8,:int16,:int32,:int64,:float32,:float64].each do |left_dtype|
    [:byte,:int8,:int16,:int32,:int64,:float32,:float64].each do |right_dtype|

      # Won't work if they're both 1-byte, due to overflow.
      next if [:byte,:int8].include?(left_dtype) && [:byte,:int8].include?(right_dtype)

      # For now, don't bother testing int-int mult.
      #next if [:int8,:int16,:int32,:int64].include?(left_dtype) && [:int8,:int16,:int32,:int64].include?(right_dtype)
      it "dense handles #{left_dtype.to_s} dot #{right_dtype.to_s} matrix multiplication" do
        #STDERR.puts "dtype=#{dtype.to_s}"
        #STDERR.puts "2"

        nary = if left_dtype.to_s =~ /complex/
                 COMPLEX_MATRIX43A_ARRAY
               else
                 MATRIX43A_ARRAY
               end

        mary = if right_dtype.to_s =~ /complex/
                 COMPLEX_MATRIX32A_ARRAY
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

  [:byte,:int8,:int16,:int32,:int64,:float32,:float64].each do |left_dtype|
    [:byte,:int8,:int16,:int32,:int64,:float32,:float64].each do |right_dtype|

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

  ALL_DTYPES.each do |dtype|
    next if integer_dtype?(dtype)
    context "#cov dtype #{dtype}" do
      before do 
        @n = NMatrix.new( [5,3], [4.0,2.0,0.60,
                                  4.2,2.1,0.59,
                                  3.9,2.0,0.58,
                                  4.3,2.1,0.62,
                                  4.1,2.2,0.63], dtype: dtype)
      end

      it "calculates variance co-variance matrix (sample)" do
        expect(@n.cov).to be_within(0.0001).of(NMatrix.new([3,3], 
          [0.025  , 0.0075, 0.00175,
           0.0075, 0.007 , 0.00135,
           0.00175, 0.00135 , 0.00043 ], dtype: dtype)
        )
      end

      it "calculates variance co-variance matrix (population)" do
        expect(@n.cov(for_sample_data: false)).to be_within(0.0001).of(NMatrix.new([3,3], 
                  [2.0000e-02, 6.0000e-03, 1.4000e-03,
                   6.0000e-03, 5.6000e-03, 1.0800e-03,
                   1.4000e-03, 1.0800e-03, 3.4400e-04], dtype: dtype)
                )
      end
    end

    context "#corr #{dtype}" do
      it "calculates the correlation matrix" do
        n = NMatrix.new([5,3], [4.0,2.0,0.60,
                                4.2,2.1,0.59,
                                3.9,2.0,0.58,
                                4.3,2.1,0.62,
                                4.1,2.2,0.63], dtype: dtype)
        expect(n.corr).to be_within(0.001).of(NMatrix.new([3,3], 
          [1.00000, 0.56695, 0.53374,
           0.56695, 1.00000, 0.77813,
           0.53374, 0.77813, 1.00000], dtype: dtype))
      end unless dtype =~ /complex/
    end

    context "#symmetric? for #{dtype} (dense)" do
      it "should return true for symmetric matrix" do 
        n = NMatrix.new([3,3], [1.00000, 0.56695, 0.53374,
                                0.56695, 1.00000, 0.77813,
                                0.53374, 0.77813, 1.00000], dtype: dtype)
        expect(n.symmetric?).to be_truthy
      end
    end
    
    context "#symmetric? for #{dtype} (list)" do
      it "should return true for symmetric matrix" do 
        n = NMatrix.new([3,3], [1.00000, 0.56695, 0.53374,
                                0.56695, 1.00000, 0.77813,
                                0.53374, 0.77813, 1.00000], dtype: dtype,
                                stype: :list)
        expect(n.symmetric?).to be_truthy
      end
      
      it "should return false for non symmetric matrix" do
        n = NMatrix.new([3,3], [1.00000, 0.56695, 1.53374,
                                0.56695, 1.00000, 0.77813,
                                0.53374, 0.77813, 1.00000], dtype: dtype,
                                stype: :list)
        expect(n.symmetric?).to be_falsy
      end
    end    

    context "#hermitian? for #{dtype} (dense)" do
      it "should return true for complex hermitian or non-complex symmetric matrix" do 
        n = NMatrix.new([3,3], [1.00000, 0.56695, 0.53374,
                                0.56695, 1.00000, 0.77813,
                                0.53374, 0.77813, 1.00000], dtype: dtype) unless dtype =~ /complex/
        n = NMatrix.new([3,3], [1.1, Complex(1.2,1.3), Complex(1.4,1.5),
                                Complex(1.2,-1.3), 1.9, Complex(1.8,1.7),
                                Complex(1.4,-1.5), Complex(1.8,-1.7), 1.3], dtype: dtype) if dtype =~ /complex/
        expect(n.hermitian?).to be_truthy
      end
    end
    
    context "#hermitian? for #{dtype} (list)" do
      it "should return true for complex hermitian or non-complex symmetric matrix" do 
        n = NMatrix.new([3,3], [1.00000, 0.56695, 0.53374,
                                0.56695, 1.00000, 0.77813,
                                0.53374, 0.77813, 1.00000], 
                                dtype: dtype, stype: :list) unless dtype =~ /complex/
        n = NMatrix.new([3,3], [1.1, Complex(1.2,1.3), Complex(1.4,1.5),
                                Complex(1.2,-1.3), 1.9, Complex(1.8,1.7),
                                Complex(1.4,-1.5), Complex(1.8,-1.7), 1.3], 
                                dtype: dtype, stype: :list) if dtype =~ /complex/
        expect(n.hermitian?).to be_truthy
      end

      it "should return false for complex non-hermitian or non-complex non-symmetric matrix" do 
        n = NMatrix.new([3,3], [1.00000, 0.56695, 0.53374,
                                0.56695, 1.00000, 1.77813,
                                0.53374, 0.77813, 1.00000], 
                                dtype: dtype, stype: :list) unless dtype =~ /complex/
        n = NMatrix.new([3,3], [1.1, Complex(1.2,1.3), Complex(1.4,-1.5),
                                Complex(1.2,-1.3), 1.9, Complex(1.8,1.7),
                                Complex(1.4,-1.5), Complex(1.8,-1.7), 1.3], 
                                dtype: dtype, stype: :list) if dtype =~ /complex/
        expect(n.hermitian?).to be_falsy
      end
    end    
    
    context "#symmetric? for #{dtype} (yale)" do
      it "should return true for symmetric matrix" do 
        n = NMatrix.new([3,3], [1.00000, 0.56695, 0.53374,
                                0.56695, 1.00000, 0.77813,
                                0.53374, 0.77813, 1.00000], dtype: dtype,
                                stype: :yale)
        expect(n.symmetric?).to be_truthy
      end
      
      it "should return false for non symmetric matrix" do
        n = NMatrix.new([3,3], [1.00000, 0.56695, 1.53374,
                                0.56695, 1.00000, 0.77813,
                                0.53374, 0.77813, 1.00000], dtype: dtype,
                                stype: :yale)
        expect(n.symmetric?).to be_falsy
      end
    end    

    context "#permute_columns for #{dtype}" do
      it "check that #permute_columns works correctly by considering every premutation of a 3x3 matrix" do
        n = NMatrix.new([3,3], [1,0,0,
                                0,2,0,
                                0,0,3], dtype: dtype)
        expect(n.permute_columns([0,1,2], {convention: :intuitive})).to eq(NMatrix.new([3,3], [1,0,0,
                                                                                              0,2,0,
                                                                                              0,0,3], dtype: dtype))
        expect(n.permute_columns([0,2,1], {convention: :intuitive})).to eq(NMatrix.new([3,3], [1,0,0,
                                                                                              0,0,2,
                                                                                              0,3,0], dtype: dtype))
        expect(n.permute_columns([1,0,2], {convention: :intuitive})).to eq(NMatrix.new([3,3], [0,1,0,
                                                                                              2,0,0,
                                                                                              0,0,3], dtype: dtype))
        expect(n.permute_columns([1,2,0], {convention: :intuitive})).to eq(NMatrix.new([3,3], [0,0,1,
                                                                                              2,0,0,
                                                                                              0,3,0], dtype: dtype))
        expect(n.permute_columns([2,0,1], {convention: :intuitive})).to eq(NMatrix.new([3,3], [0,1,0,
                                                                                              0,0,2,
                                                                                              3,0,0], dtype: dtype))
        expect(n.permute_columns([2,1,0], {convention: :intuitive})).to eq(NMatrix.new([3,3], [0,0,1,
                                                                                              0,2,0,
                                                                                              3,0,0], dtype: dtype))
        expect(n.permute_columns([0,1,2], {convention: :lapack})).to eq(NMatrix.new([3,3], [1,0,0,
                                                                                           0,2,0,
                                                                                           0,0,3], dtype: dtype))
        expect(n.permute_columns([0,2,2], {convention: :lapack})).to eq(NMatrix.new([3,3], [1,0,0,
                                                                                           0,0,2,
                                                                                           0,3,0], dtype: dtype))
        expect(n.permute_columns([1,1,2], {convention: :lapack})).to eq(NMatrix.new([3,3], [0,1,0,
                                                                                           2,0,0,
                                                                                           0,0,3], dtype: dtype))
        expect(n.permute_columns([1,2,2], {convention: :lapack})).to eq(NMatrix.new([3,3], [0,0,1,
                                                                                           2,0,0,
                                                                                           0,3,0], dtype: dtype))
        expect(n.permute_columns([2,2,2], {convention: :lapack})).to eq(NMatrix.new([3,3], [0,1,0,
                                                                                           0,0,2,
                                                                                           3,0,0], dtype: dtype))
        expect(n.permute_columns([2,1,2], {convention: :lapack})).to eq(NMatrix.new([3,3], [0,0,1,
                                                                                           0,2,0,
                                                                                           3,0,0], dtype: dtype))
      end
      it "additional tests for  #permute_columns with convention :intuitive" do
        m = NMatrix.new([1,4], [0,1,2,3], dtype: dtype)
        perm = [1,0,3,2]
        expect(m.permute_columns(perm, {convention: :intuitive})).to eq(NMatrix.new([1,4], perm, dtype: dtype))

        m = NMatrix.new([1,5], [0,1,2,3,4], dtype: dtype)
        perm = [1,0,4,3,2]
        expect(m.permute_columns(perm, {convention: :intuitive})).to eq(NMatrix.new([1,5], perm, dtype: dtype))

        m = NMatrix.new([1,6], [0,1,2,3,4,5], dtype: dtype)
        perm = [2,4,1,0,5,3]
        expect(m.permute_columns(perm, {convention: :intuitive})).to eq(NMatrix.new([1,6], perm, dtype: dtype))

        m = NMatrix.new([1,7], [0,1,2,3,4,5,6], dtype: dtype)
        perm = [1,3,5,6,0,2,4]
        expect(m.permute_columns(perm, {convention: :intuitive})).to eq(NMatrix.new([1,7], perm, dtype: dtype))

        m = NMatrix.new([1,8], [0,1,2,3,4,5,6,7], dtype: dtype)
        perm = [6,7,5,4,1,3,0,2]
        expect(m.permute_columns(perm, {convention: :intuitive})).to eq(NMatrix.new([1,8], perm, dtype: dtype))
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

      it "solves linear equation for #{dtype} (non-symmetric matrix)" do
        a = NMatrix.new [3,3], [1,1,1, -1,0,1, 3,4,6], dtype: dtype
        b = NMatrix.new [3,1], [6,2,29], dtype: dtype

        err = case dtype
                when :float32, :complex64
                  1e-5
                else
                  1e-14
              end

        expect(a.solve(b)).to be_within(err).of(NMatrix.new([3,1], [1,2,3], dtype: dtype))
      end

      it "solves linear equation for dtype #{dtype} (non-vector rhs)" do
        a = NMatrix.new [3,3], [1,0,0, -1,0,1, 2,1,1], dtype: dtype
        b = NMatrix.new [3,2], [1,0, 1,2, 4,2], dtype: dtype

        expect(a.solve(b)).to eq(NMatrix.new [3,2], [1,0, 0,0, 2,2], dtype: dtype)
      end
    end

    FLOAT_DTYPES.each do |dtype|
      context "when form: :lower_tri" do
        let(:a) { NMatrix.new([3,3], [1, 0, 0, 2, 0.5, 0, 3, 3, 9], dtype: dtype) }

        it "solves a lower triangular linear system A * x = b with vector b" do
          b = NMatrix.new([3,1], [1,2,3], dtype: dtype)
          x = a.solve(b, form: :lower_tri)
          r = a.dot(x) - b
          expect(r.abs.max).to be_within(1e-6).of(0.0)
        end

        it "solves a lower triangular linear system A * X = B with narrow B" do
          b = NMatrix.new([3,2], [1,2,3,4,5,6], dtype: dtype)
          x = a.solve(b, form: :lower_tri)
          r = (a.dot(x) - b).abs.to_flat_a
          expect(r.max).to be_within(1e-6).of(0.0)
        end

        it "solves a lower triangular linear system A * X = B with wide B" do
          b = NMatrix.new([3,5], (1..15).to_a, dtype: dtype)
          x = a.solve(b, form: :lower_tri)
          r = (a.dot(x) - b).abs.to_flat_a
          expect(r.max).to be_within(1e-6).of(0.0)
        end
      end

      context "when form: :upper_tri" do
        let(:a) { NMatrix.new([3,3], [3, 2, 1, 0, 2, 0.5, 0, 0, 9], dtype: dtype) }

        it "solves an upper triangular linear system A * x = b with vector b" do
          b = NMatrix.new([3,1], [1,2,3], dtype: dtype)
          x = a.solve(b, form: :upper_tri)
          r = a.dot(x) - b
          expect(r.abs.max).to be_within(1e-6).of(0.0)
        end

        it "solves an upper triangular linear system A * X = B with narrow B" do
          b = NMatrix.new([3,2], [1,2,3,4,5,6], dtype: dtype)
          x = a.solve(b, form: :upper_tri)
          r = (a.dot(x) - b).abs.to_flat_a
          expect(r.max).to be_within(1e-6).of(0.0)
        end

        it "solves an upper triangular linear system A * X = B with a wide B" do
          b = NMatrix.new([3,5], (1..15).to_a, dtype: dtype)
          x = a.solve(b, form: :upper_tri)
          r = (a.dot(x) - b).abs.to_flat_a
          expect(r.max).to be_within(1e-6).of(0.0)
        end
      end

      context "when form: :pos_def" do
        let(:a) { NMatrix.new([3,3], [4, 1, 2, 1, 5, 3, 2, 3, 6], dtype: dtype) }

        it "solves a linear system A * X = b with positive definite A and vector b" do
          b = NMatrix.new([3,1], [6,4,8], dtype: dtype)
          begin
            x = a.solve(b, form: :pos_def)
            expect(x).to be_within(1e-6).of(NMatrix.new([3,1], [1,0,1], dtype: dtype))
          rescue NotImplementedError
            "Suppressing a NotImplementedError when the lapacke or atlas plugin is not available"
          end
        end
      
        it "solves a linear system A * X = B with positive definite A and matrix B" do
          b = NMatrix.new([3,2], [8,3,14,13,14,19], dtype: dtype)
          begin
            x = a.solve(b, form: :pos_def)
            expect(x).to be_within(1e-6).of(NMatrix.new([3,2], [1,-1,2,1,1,3], dtype: dtype))
          rescue NotImplementedError
            "Suppressing a NotImplementedError when the lapacke or atlas plugin is not available"
          end
        end
      end
    end
  end

  context "#hessenberg" do
    FLOAT_DTYPES.each do |dtype|
      context dtype do
        before do
          @n = NMatrix.new [5,5], 
            [0, 2, 0, 1, 1,
             2, 2, 3, 2, 2,
             4,-3, 0, 1, 3,
             6, 1,-6,-5, 4,
             5, 6, 4, 1, 5], dtype: dtype
        end

        it "transforms a matrix to Hessenberg form" do
          expect(@n.hessenberg).to be_within(0.0001).of(NMatrix.new([5,5],    
            [0.00000,-1.66667, 0.79432,-0.45191,-1.54501,
            -9.00000, 2.95062,-6.89312, 3.22250,-0.19012,
             0.00000,-8.21682,-0.57379, 5.26966,-1.69976,
             0.00000, 0.00000,-3.74630,-0.80893, 3.99708,
             0.00000, 0.00000, 0.00000, 0.04102, 0.43211], dtype: dtype))
        end
      end
    end
  end

  ALL_DTYPES.each do |dtype|
    [:dense, :yale].each do |stype|
      answer_dtype = integer_dtype?(dtype) ? :int64 : dtype
      next if dtype == :byte
      
      context "#pow #{dtype} #{stype}" do
        before do 
          @n = NMatrix.new [4,4], [0, 2, 0, 1,
                                  2, 2, 3, 2,
                                  4,-3, 0, 1,
                                  6, 1,-6,-5], dtype: dtype, stype: stype
        end

        it "raises a square matrix to even power" do
          expect(@n.pow(4)).to eq(NMatrix.new([4,4], [292, 28,-63, -42, 
                                                     360, 96, 51, -14,
                                                     448,-231,-24,-87,
                                                   -1168, 595,234, 523], 
                                                   dtype: answer_dtype, 
                                                   stype: stype))
        end

        it "raises a square matrix to odd power" do
          expect(@n.pow(9)).to eq(NMatrix.new([4,4],[-275128,  279917, 176127, 237451,
                                                    -260104,  394759, 166893,  296081,
                                                    -704824,  285700, 186411,  262002,
                                                    3209256,-1070870,-918741,-1318584],
                                                    dtype: answer_dtype, stype: stype))
        end

        it "raises a sqaure matrix to negative power" do
          expect(@n.pow(-3)).to be_within(0.00001).of (NMatrix.new([4,4],
            [1.0647e-02, 4.2239e-04,-6.2281e-05, 2.7680e-03,
            -1.6415e-02, 2.1296e-02, 1.0718e-02, 4.8589e-03,   
             8.6956e-03,-8.6569e-03, 2.8993e-02, 7.2015e-03,
             5.0034e-02,-1.7500e-02,-3.6777e-02,-1.2128e-02], dtype: answer_dtype, 
             stype: stype)) 
        end unless stype =~ /yale/ or dtype == :object or ALL_DTYPES.grep(/int/).include? dtype

        it "raises a square matrix to zero" do
          expect(@n.pow(0)).to eq(NMatrix.eye([4,4], dtype: answer_dtype, 
            stype: stype))
        end

        it "raises a square matrix to one" do
          expect(@n.pow(1)).to eq(@n)
        end
      end
    end
  end

  ALL_DTYPES.each do |dtype|
    [:dense, :yale].each do |stype|
      context "#kron_prod #{dtype} #{stype}" do
        before do 
          @a = NMatrix.new([2,2], [1,2,
                                   3,4], dtype: dtype, stype: stype)
          @b = NMatrix.new([2,3], [1,1,1,
                                   1,1,1], dtype: dtype, stype: stype)
          @c = NMatrix.new([4,6], [1, 1, 1, 2, 2, 2,
                                   1, 1, 1, 2, 2, 2,
                                   3, 3, 3, 4, 4, 4,
                                   3, 3, 3, 4, 4, 4], dtype: dtype, stype: stype) 
        end
        it "Compute the Kronecker product of two NMatrix" do
          expect(@a.kron_prod(@b)).to eq(@c)
        end
      end
    end
  end

  context "determinants" do
    ALL_DTYPES.each do |dtype|
      next if dtype == :object
      context dtype do
        before do
          @a = NMatrix.new([2,2], [1,2,
                                   3,4], dtype: dtype)
          @b = NMatrix.new([3,3], [1,2,3,
                                   5,0,1,
                                   4,1,3], dtype: dtype)
          @c = NMatrix.new([4,4], [1, 0, 1, 1,
                                   1, 2, 3, 1,
                                   3, 3, 3, 1,
                                   1, 2, 3, 4], dtype: dtype)
          @err = case dtype
                  when :float32, :complex64
                    1e-6
                  when :float64, :complex128
                    1e-14
                  else
                    1e-64 # FIXME: should be 0, but be_within(0) does not work.
                end
        end
        it "computes the determinant of 2x2 matrix" do
          expect(@a.det).to be_within(@err).of(-2)
        end
        it "computes the determinant of 3x3 matrix" do
          expect(@b.det).to be_within(@err).of(-8)
        end
        it "computes the determinant of 4x4 matrix" do
          expect(@c.det).to be_within(@err).of(-18)
        end
        it "computes the exact determinant of 2x2 matrix" do
          if dtype == :byte
            expect{@a.det_exact}.to raise_error(DataTypeError)
          else
            expect(@a.det_exact).to be_within(@err).of(-2)
          end
        end
        it "computes the exact determinant of 3x3 matrix" do
          if dtype == :byte
            expect{@a.det_exact}.to raise_error(DataTypeError)
          else
            expect(@b.det_exact).to be_within(@err).of(-8)
          end
        end
      end
    end
  end
end
