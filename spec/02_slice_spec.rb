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
require_relative 'spec_helper'

describe "Slice operation" do
  include RSpec::Longrun::DSL

  [:dense, :list, :yale].each do |stype|
    context "for #{stype}" do
        #GC.start # don't have to do this, but it helps to make sure we've cleaned up our pointers properly.
      let(:stype_matrix) { create_matrix(stype) }

      it "should correctly return a row of a reference-slice" do
        n = create_rectangular_matrix(stype)
        stype_matrix = n[1..4,1..3]
        expect(stype_matrix.row(1, :copy)).to eq(stype_matrix.row(1, :reference))
        expect(stype_matrix.row(1, :copy).to_flat_array).to eq([12,13,0])
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
          stype_matrix[2,1..2].each_stored_with_indices do |v,i,j|
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
          stype_matrix[2,0..2].each_stored_with_indices do |v,i,j|
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
          stype_matrix[0..1,1..2].each_stored_with_indices do |v,i,j|
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
        a = stype_matrix[0..1, 0..1]
        b = stype_matrix.slice(0..1, 0..1)
        expect(stype_matrix.is_ref?).to be false
        expect(a.is_ref?).to be true
        expect(b.is_ref?).to be false
      end

      it "reference should compare with non-reference" do
        expect(stype_matrix.slice(1..2,0..1)).to eq(stype_matrix[1..2, 0..1])
        expect(stype_matrix[1..2,0..1]).to eq(stype_matrix.slice(1..2, 0..1))
        expect(stype_matrix[1..2,0..1]).to eq(stype_matrix[1..2, 0..1])
      end

      context "with copying" do
        it 'should return an NMatrix' do
          n = stype_matrix.slice(0..1,0..1)
          expect(nm_eql(n, NMatrix.new([2,2], [0,1,3,4], dtype: :int32))).to be true
        end

        it 'should return a copy of 2x2 matrix to self elements' do
          n = stype_matrix.slice(1..2,0..1)
          expect(n.shape).to eql([2,2])

          expect(n[1,1]).to eq(stype_matrix[2,1])
          n[1,1] = -9
          expect(stype_matrix[2,1]).to eql(7)
        end

        it 'should return a 1x2 matrix without refs to self elements' do
          n = stype_matrix.slice(0,1..2)
          expect(n.shape).to eql([1,2])

          expect(n[0]).to eq(stype_matrix[0,1])
          expect(n[1]).to eq(stype_matrix[0,2])
          n[0] = -9
          expect(stype_matrix[0,1]).to eql(1)
          expect(stype_matrix[0,2]).to eql(2)
        end

        it 'should return a 2x1 matrix without refs to self elements' do
          stype_matrix.extend NMatrix::YaleFunctions

          n = stype_matrix.slice(0..1,1)
          expect(n.shape).to eql([2,1])

          expect(n[0]).to eq(stype_matrix[0,1])
          expect(n[1]).to eq(stype_matrix[1,1])
          n[0] = -9
          expect(stype_matrix[0,1]).to eql(1)
          expect(stype_matrix[1,1]).to eql(4)
        end

        it 'should be correct slice for range 0..2 and 0...3' do
          expect(stype_matrix.slice(0..2,0..2)).to eq(stype_matrix.slice(0...3,0...3))
        end

        [:dense, :list, :yale].each do |cast_type|
          it "should cast copied slice from #{stype.upcase} to #{cast_type.upcase}" do
            expect(nm_eql(stype_matrix.slice(1..2, 1..2).cast(cast_type, :int32), stype_matrix.slice(1..2,1..2))).to be true
            expect(nm_eql(stype_matrix.slice(0..1, 1..2).cast(cast_type, :int32), stype_matrix.slice(0..1,1..2))).to be true
            expect(nm_eql(stype_matrix.slice(1..2, 0..1).cast(cast_type, :int32), stype_matrix.slice(1..2,0..1))).to be true
            expect(nm_eql(stype_matrix.slice(0..1, 0..1).cast(cast_type, :int32), stype_matrix.slice(0..1,0..1))).to be true

            # Non square
            expect(nm_eql(stype_matrix.slice(0..2, 1..2).cast(cast_type, :int32), stype_matrix.slice(0..2,1..2))).to be true
            #require 'pry'
            #binding.pry if cast_type == :yale
            expect(nm_eql(stype_matrix.slice(1..2, 0..2).cast(cast_type, :int32), stype_matrix.slice(1..2,0..2))).to be true

            # Full
            expect(nm_eql(stype_matrix.slice(0..2, 0..2).cast(cast_type, :int32), stype_matrix)).to be true
          end
        end
      end

      # Yale:
      #context "by copy" do
        #it "should correctly preserve zeros" do
        #  stype_matrix = NMatrix.new(:yale, 3, :int64)
        #  column_slice = stype_matrix.column(2, :copy)
        #  column_slice[0].should == 0
        #  column_slice[1].should == 0
        #  column_slice[2].should == 0
        #end
      #end

      context "by reference" do
        it 'should return an NMatrix' do
          n = stype_matrix[0..1,0..1]
          expect(nm_eql(n, NMatrix.new([2,2], [0,1,3,4], dtype: :int32))).to be true
        end

        it 'should return a 2x2 matrix with refs to self elements' do
          n = stype_matrix[1..2,0..1]
          expect(n.shape).to eql([2,2])

          expect(n[0,0]).to eq(stype_matrix[1,0])
          n[0,0] = -9
          expect(stype_matrix[1,0]).to eql(-9)
        end

        it 'should return a 1x2 vector with refs to self elements' do
          n = stype_matrix[0,1..2]
          expect(n.shape).to eql([1,2])

          expect(n[0]).to eq(stype_matrix[0,1])
          n[0] = -9
          expect(stype_matrix[0,1]).to eql(-9)
        end

        it 'should return a 2x1 vector with refs to self elements' do
          n = stype_matrix[0..1,1]
          expect(n.shape).to eql([2,1])

          expect(n[0]).to eq(stype_matrix[0,1])
          n[0] = -9
          expect(stype_matrix[0,1]).to eql(-9)
        end

        it 'should slice again' do
          n = stype_matrix[1..2, 1..2]
          expect(nm_eql(n[1,0..1], NVector.new(2, [7,8], dtype: :int32).transpose)).to be true
        end

        it 'should be correct slice for range 0..2 and 0...3' do
          expect(stype_matrix[0..2,0..2]).to eq(stype_matrix[0...3,0...3])
        end

        it 'should correctly handle :* slice notation' do
          expect(stype_matrix[:*,0]).to eq stype_matrix[0...stype_matrix.shape[0], 0]
        end

        if stype == :dense
          [:byte,:int8,:int16,:int32,:int64,:float32,:float64].each do |left_dtype|
            [:byte,:int8,:int16,:int32,:int64,:float32,:float64].each do |right_dtype|

              # Won't work if they're both 1-byte, due to overflow.
              next if [:byte,:int8].include?(left_dtype) && [:byte,:int8].include?(right_dtype)

              # For now, don't bother testing int-int mult.
              #next if [:int8,:int16,:int32,:int64].include?(left_dtype) && [:int8,:int16,:int32,:int64].include?(right_dtype)
              it "handles #{left_dtype.to_s} dot #{right_dtype.to_s} matrix multiplication" do
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
              expect(stype_matrix[0...3,0].transpose).to eq NMatrix[[0, 3, 6]]
              expect(stype_matrix[0...3,1].transpose).to eq NMatrix[[1, 4, 7]]
              expect(stype_matrix[0...3,2].transpose).to eq NMatrix[[2, 5, 8]]
              expect(stype_matrix[0,0...3].transpose).to eq NMatrix[[0], [1], [2]]
              expect(stype_matrix[1,0...3].transpose).to eq NMatrix[[3], [4], [5]]
              expect(stype_matrix[2,0...3].transpose).to eq NMatrix[[6], [7], [8]]
              expect(stype_matrix[1..2,1..2].transpose).to eq NMatrix[[4, 7], [5, 8]]
            end

            it "adds slices" do
              expect(NMatrix[[0,0,0]] + stype_matrix[1,0..2]).to eq NMatrix[[3, 4, 5]]
            end

            it "scalar adds to slices" do
              expect(stype_matrix[1,0..2]+1).to eq NMatrix[[4, 5, 6]]
            end

            it "compares slices to scalars" do
              (stype_matrix[1, 0..2] > 2).each { |e| expect(e != 0).to be true }
            end

            it "iterates only over elements in the slice" do
              els = []
              stype_matrix[1, 0..2].each { |e| els << e }
              expect(els.size).to eq 3
              expect(els[0]).to eq 3
              expect(els[1]).to eq 4
              expect(els[2]).to eq 5
            end

            it "iterates with index only over elements in the slice" do
              els = []
              stype_matrix[1, 0..2].each_stored_with_indices { |a| els << a }
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
              n = stype_matrix[1..2,0..1]
            end
            GC.start
          end

          step "reference slice of casted-copy" do
            expect(stype_matrix).to eq(NMatrix.new([3,3], (0..9).to_a, dtype: :int32).cast(stype, :int32))
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
            expect(nm_eql(stype_matrix[1..2, 1..2].cast(cast_type), stype_matrix[1..2,1..2])).to be true
            expect(nm_eql(stype_matrix[0..1, 1..2].cast(cast_type), stype_matrix[0..1,1..2])).to be true
            expect(nm_eql(stype_matrix[1..2, 0..1].cast(cast_type), stype_matrix[1..2,0..1])).to be true
            expect(nm_eql(stype_matrix[0..1, 0..1].cast(cast_type), stype_matrix[0..1,0..1])).to be true
          end

          it "should cast a rectangular reference-slice from #{stype.upcase} to #{cast_type.upcase}" do
            # Non square
            expect(nm_eql(stype_matrix[0..2, 1..2].cast(cast_type), stype_matrix[0..2,1..2])).to be true # FIXME: memory problem.
            expect(nm_eql(stype_matrix[1..2, 0..2].cast(cast_type), stype_matrix[1..2,0..2])).to be true # this one is fine
          end

          it "should cast a square full-matrix reference-slice from #{stype.upcase} to #{cast_type.upcase}" do
            expect(nm_eql(stype_matrix[0..2, 0..2].cast(cast_type), stype_matrix)).to be true
          end
        end
      end
    end
  end
end
