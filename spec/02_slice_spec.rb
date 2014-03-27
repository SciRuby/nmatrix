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
# == 02_slice_spec.rb
#
# Test of slice operations. High priority tests since reference
# slicing is needed for pretty_print.
#
require 'pry'
require File.dirname(__FILE__) + "/spec_helper.rb"

describe "Slice operation" do
  include RSpec::Longrun::DSL

  [:dense, :list, :yale].each do |stype|
    context "for #{stype}" do
      before :each do
        #GC.start # don't have to do this, but it helps to make sure we've cleaned up our pointers properly.
        @m = create_matrix(stype)
      end

      #after :each do
      #  GC.start
      #end

      it "should correctly return a row of a reference-slice" do
        @n = create_rectangular_matrix(stype)
        @m = @n[1..4,1..3]
        expect(@m.row(1, :copy)).to eq(@m.row(1, :reference))
        expect(@m.row(1, :copy).to_flat_array).to eq([12,13,0])
      end

      if stype == :yale
        it "should binary search for the left boundary of a partial row of stored indices correctly" do
          n = NMatrix.new(10, stype: :yale, dtype: :int32)
          n[3,0] = 1
          #n[3,2] = 2
          n[3,3] = 3
          n[3,4] = 4
          n[3,6] = 5
          n[3,8] = 6
          n[3,9] = 7
          vs = []
          is = []
          js = []

          n[3,1..9].each_stored_with_indices do |v,i,j|
            vs << v
            is << i
            js << j
          end

          expect(vs).to eq([3,4,5,6,7])
          expect(js).to eq([2,3,5,7,8])
          expect(is).to eq([0,0,0,0,0])
        end
      elsif stype == :list
        it "should iterate across a partial row of stored indices" do
          vs = []
          is = []
          js = []

          STDERR.puts("now") if stype == :yale
          @m[2,1..2].each_stored_with_indices do |v,i,j|
            vs << v
            is << i
            js << j
          end

          expect(vs).to eq([7,8])
          expect(is).to eq([0,0])
          expect(js).to eq([0,1])
        end
      end

      unless stype == :dense
        it "should iterate across a row of stored indices" do

          vs = []
          is = []
          js = []
          @m[2,0..2].each_stored_with_indices do |v,i,j|
            vs << v
            is << i
            js << j
          end
          expect(vs).to eq(stype == :yale ? [8,6,7] : [6,7,8])
          expect(is).to eq([0,0,0])
          expect(js).to eq(stype == :yale ? [2,0,1] : [0,1,2])
        end

        it "should iterate across a submatrix of stored indices" do
          vs = []
          is = []
          js = []
          @m[0..1,1..2].each_stored_with_indices do |v,i,j|
            vs << v
            is << i
            js << j
          end

          expect(vs).to eq(stype == :yale ? [4,1,2,5] : [1,2,4,5])
          expect(is).to eq(stype == :yale ? [1,0,0,1] : [0,0,1,1])
          expect(js).to eq(stype == :yale ? [0,0,1,1] : [0,1,0,1])
        end
      end

      it "should return correct supershape" do
        x = NMatrix.random([10,12])
        y = x[0...8,5...12]
        expect(y.shape).to eq([8,7])
        expect(y.supershape).to eq([10,12])
      end

      it "should have #is_ref? method" do
        a = @m[0..1, 0..1]
        b = @m.slice(0..1, 0..1)
        expect(@m.is_ref?).to be_false
        expect(a.is_ref?).to be_true
        expect(b.is_ref?).to be_false
      end

      it "reference should compare with non-reference" do
        expect(@m.slice(1..2,0..1)).to eq(@m[1..2, 0..1])
        expect(@m[1..2,0..1]).to eq(@m.slice(1..2, 0..1))
        expect(@m[1..2,0..1]).to eq(@m[1..2, 0..1])
      end

      context "with copying" do
        it 'should return an NMatrix' do
          n = @m.slice(0..1,0..1)
          expect(nm_eql(n, NMatrix.new([2,2], [0,1,3,4], dtype: :int32))).to be_true
        end

        it 'should return a copy of 2x2 matrix to self elements' do
          n = @m.slice(1..2,0..1)
          expect(n.shape).to eql([2,2])

          expect(n[1,1]).to eq(@m[2,1])
          n[1,1] = -9
          expect(@m[2,1]).to eql(7)
        end

        it 'should return a 1x2 matrix without refs to self elements' do
          n = @m.slice(0,1..2)
          expect(n.shape).to eql([1,2])

          expect(n[0]).to eq(@m[0,1])
          expect(n[1]).to eq(@m[0,2])
          n[0] = -9
          expect(@m[0,1]).to eql(1)
          expect(@m[0,2]).to eql(2)
        end

        it 'should return a 2x1 matrix without refs to self elements' do
          @m.extend NMatrix::YaleFunctions

          n = @m.slice(0..1,1)
          expect(n.shape).to eql([2,1])

          expect(n[0]).to eq(@m[0,1])
          expect(n[1]).to eq(@m[1,1])
          n[0] = -9
          expect(@m[0,1]).to eql(1)
          expect(@m[1,1]).to eql(4)
        end

        it 'should be correct slice for range 0..2 and 0...3' do
          expect(@m.slice(0..2,0..2)).to eq(@m.slice(0...3,0...3))
        end

        [:dense, :list, :yale].each do |cast_type|
          it "should cast copied slice from #{stype.upcase} to #{cast_type.upcase}" do
            expect(nm_eql(@m.slice(1..2, 1..2).cast(cast_type, :int32), @m.slice(1..2,1..2))).to be_true
            expect(nm_eql(@m.slice(0..1, 1..2).cast(cast_type, :int32), @m.slice(0..1,1..2))).to be_true
            expect(nm_eql(@m.slice(1..2, 0..1).cast(cast_type, :int32), @m.slice(1..2,0..1))).to be_true
            expect(nm_eql(@m.slice(0..1, 0..1).cast(cast_type, :int32), @m.slice(0..1,0..1))).to be_true

            # Non square
            expect(nm_eql(@m.slice(0..2, 1..2).cast(cast_type, :int32), @m.slice(0..2,1..2))).to be_true
            #require 'pry'
            #binding.pry if cast_type == :yale
            expect(nm_eql(@m.slice(1..2, 0..2).cast(cast_type, :int32), @m.slice(1..2,0..2))).to be_true

            # Full
            expect(nm_eql(@m.slice(0..2, 0..2).cast(cast_type, :int32), @m)).to be_true
          end
        end
      end

      # Yale:
      #context "by copy" do
        #it "should correctly preserve zeros" do
        #  @m = NMatrix.new(:yale, 3, :int64)
        #  column_slice = @m.column(2, :copy)
        #  column_slice[0].should == 0
        #  column_slice[1].should == 0
        #  column_slice[2].should == 0
        #end
      #end

      context "by reference" do
        it 'should return an NMatrix' do
          n = @m[0..1,0..1]
          expect(nm_eql(n, NMatrix.new([2,2], [0,1,3,4], dtype: :int32))).to be_true
        end

        it 'should return a 2x2 matrix with refs to self elements' do
          n = @m[1..2,0..1]
          expect(n.shape).to eql([2,2])

          expect(n[0,0]).to eq(@m[1,0])
          n[0,0] = -9
          expect(@m[1,0]).to eql(-9)
        end

        it 'should return a 1x2 vector with refs to self elements' do
          n = @m[0,1..2]
          expect(n.shape).to eql([1,2])

          expect(n[0]).to eq(@m[0,1])
          n[0] = -9
          expect(@m[0,1]).to eql(-9)
        end

        it 'should return a 2x1 vector with refs to self elements' do
          n = @m[0..1,1]
          expect(n.shape).to eql([2,1])

          expect(n[0]).to eq(@m[0,1])
          n[0] = -9
          expect(@m[0,1]).to eql(-9)
        end

        it 'should slice again' do
          n = @m[1..2, 1..2]
          expect(nm_eql(n[1,0..1], NVector.new(2, [7,8], dtype: :int32).transpose)).to be_true
        end

        it 'should be correct slice for range 0..2 and 0...3' do
          expect(@m[0..2,0..2]).to eq(@m[0...3,0...3])
        end

        it 'should correctly handle :* slice notation' do
          expect(@m[:*,0]).to eq @m[0...@m.shape[0], 0]
        end

        if stype == :dense
          [:byte,:int8,:int16,:int32,:int64,:float32,:float64,:rational64,:rational128].each do |left_dtype|
            [:byte,:int8,:int16,:int32,:int64,:float32,:float64,:rational64,:rational128].each do |right_dtype|

              # Won't work if they're both 1-byte, due to overflow.
              next if [:byte,:int8].include?(left_dtype) && [:byte,:int8].include?(right_dtype)

              # For now, don't bother testing int-int mult.
              #next if [:int8,:int16,:int32,:int64].include?(left_dtype) && [:int8,:int16,:int32,:int64].include?(right_dtype)
              it "handles #{left_dtype.to_s} dot #{right_dtype.to_s} matrix multiplication" do
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

                n = NMatrix.new([4,3], nary, dtype: left_dtype)[1..3,1..2]
                m = NMatrix.new([3,2], mary, dtype: right_dtype)[1..2,0..1]

                r = n.dot m
                expect(r.shape).to eql([3,2])

                expect(r[0,0]).to eq(219.0)
                expect(r[0,1]).to eq(185.0)
                expect(r[1,0]).to eq(244.0)
                expect(r[1,1]).to eq(205.0)
                expect(r[2,0]).to eq(42.0)
                expect(r[2,1]).to eq(35.0)

              end
            end
          end

          context "operations" do

            it "correctly transposes slices" do
              expect(@m[0...3,0].transpose).to eq NMatrix[[0, 3, 6]]
              expect(@m[0...3,1].transpose).to eq NMatrix[[1, 4, 7]]
              expect(@m[0...3,2].transpose).to eq NMatrix[[2, 5, 8]]
              expect(@m[0,0...3].transpose).to eq NMatrix[[0], [1], [2]]
              expect(@m[1,0...3].transpose).to eq NMatrix[[3], [4], [5]]
              expect(@m[2,0...3].transpose).to eq NMatrix[[6], [7], [8]]
            end

            it "adds slices" do
              expect(NMatrix[[0,0,0]] + @m[1,0..2]).to eq NMatrix[[3, 4, 5]]
            end

            it "scalar adds to slices" do
              expect(@m[1,0..2]+1).to eq NMatrix[[4, 5, 6]]
            end

            it "compares slices to scalars" do
              (@m[1, 0..2] > 2).each { |e| expect(e != 0).to be_true }
            end

            it "iterates only over elements in the slice" do
              els = []
              @m[1, 0..2].each { |e| els << e }
              expect(els.size).to eq 3
              expect(els[0]).to eq 3
              expect(els[1]).to eq 4
              expect(els[2]).to eq 5
            end

            it "iterates with index only over elements in the slice" do
              els = []
              @m[1, 0..2].each_stored_with_indices { |a| els << a }
              expect(els.size).to eq 3
              expect(els[0]).to eq [3, 0, 0]
              expect(els[1]).to eq [4, 0, 1]
              expect(els[2]).to eq [5, 0, 2]
            end

          end

        end

        example 'should be cleaned up by garbage collector without errors'  do
          step "reference slice" do
            1.times do
              n = @m[1..2,0..1]
            end
            GC.start
          end

          step "reference slice of casted-copy" do
            expect(@m).to eq(NMatrix.new([3,3], (0..9).to_a, dtype: :int32).cast(stype, :int32))
            n = nil
            1.times do
              m = NMatrix.new([2,2], [1,2,3,4]).cast(stype, :int32)
              n = m[0..1,0..1]
            end
            GC.start
            expect(n).to eq(NMatrix.new([2,2], [1,2,3,4]).cast(stype, :int32))
          end
        end

        [:dense, :list, :yale].each do |cast_type|
          it "should cast a square reference-slice from #{stype.upcase} to #{cast_type.upcase}" do
            expect(nm_eql(@m[1..2, 1..2].cast(cast_type), @m[1..2,1..2])).to be_true
            expect(nm_eql(@m[0..1, 1..2].cast(cast_type), @m[0..1,1..2])).to be_true
            expect(nm_eql(@m[1..2, 0..1].cast(cast_type), @m[1..2,0..1])).to be_true
            expect(nm_eql(@m[0..1, 0..1].cast(cast_type), @m[0..1,0..1])).to be_true
          end

          it "should cast a rectangular reference-slice from #{stype.upcase} to #{cast_type.upcase}" do
            # Non square
            expect(nm_eql(@m[0..2, 1..2].cast(cast_type), @m[0..2,1..2])).to be_true # FIXME: memory problem.
            expect(nm_eql(@m[1..2, 0..2].cast(cast_type), @m[1..2,0..2])).to be_true # this one is fine
          end

          it "should cast a square full-matrix reference-slice from #{stype.upcase} to #{cast_type.upcase}" do
            expect(nm_eql(@m[0..2, 0..2].cast(cast_type), @m)).to be_true
          end
        end
      end

    end
  end
end
