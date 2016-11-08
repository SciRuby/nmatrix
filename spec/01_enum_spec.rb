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
# == 01_enum_spec.rb
#
# Enumerator tests for NMatrix. These should load early, as they
# test functionality essential to matrix printing.
#
require_relative 'spec_helper'

describe "NMatrix enumeration for" do
  [:dense, :yale, :list].each do |stype|
    context stype do
      let(:n) { create_rectangular_matrix(stype) }
      let(:m) { n[1..4,1..3] }

      if stype == :yale
        it "should iterate properly along each row of a slice" do
          vv = []
          ii = []
          jj = []
          m.extend NMatrix::YaleFunctions
          m.each_row do |row|
            row.each_with_indices do |v,i,j|
              vv << v
              ii << i
              jj << j
            end
          end

          expect(vv).to eq([7,8,9, 12,13,0, 0,0,0, 0,17,18])
          expect(ii).to eq([0]*12)
          expect(jj).to eq([0,1,2]*4)
        end

        it "should iterate along diagonal portion of A array" do
          vv = []
          ii = []
          jj = []
          n.send :__yale_stored_diagonal_each_with_indices__ do |v,i,j|
            vv << v
            ii << i
            jj << j
          end
          expect(vv).to eq([1,7,13,0,19])
          expect(ii).to eq([0,1,2,3,4])
          expect(jj).to eq(ii)
        end

        it "should iterate along non-diagonal portion of A array" do
          vv = []
          ii = []
          jj = []
          n.send :__yale_stored_nondiagonal_each_with_indices__ do |v,i,j|
            vv << v
            ii << i
            jj << j
          end

          expect(vv).to eq([2,3,4,5,  6,8,9,10,  11,12,14,15,  16,17,18,20])
          expect(ii).to eq([[0]*4, [1]*4, [2]*4, [4]*4].flatten)
          expect(jj).to eq([1,2,3,4,  0,2,3,5,   0,1,4,5,      0,2,3,5])
        end

        it "should iterate along a sliced diagonal portion of an A array" do
          m = n[0..3,1..3]
          vv = []
          ii = []
          jj = []
          m.send :__yale_stored_diagonal_each_with_indices__ do |v,i,j|
            vv << v
            ii << i
            jj << j
          end
          expect(vv).to eq([7,13,0])
          expect(ii).to eq([1,2,3])
          expect(jj).to eq([0,1,2])
        end

        it "should iterate along a sliced non-diagonal portion of a sliced A array" do
          vv = []
          ii = []
          jj = []
          n.extend NMatrix::YaleFunctions
          m.extend NMatrix::YaleFunctions
          m.send :__yale_stored_nondiagonal_each_with_indices__ do |v,i,j|
            vv << v
            ii << i
            jj << j
          end

          expect(ii).to eq([0,0, 1,   3,3 ])
          expect(jj).to eq([1,2, 0,   1,2 ])
          expect(vv).to eq([8,9, 12, 17,18])
        end

        it "should visit each stored element of the matrix in order by indices" do
          vv = []
          ii = []
          jj = []
          n.each_ordered_stored_with_indices do |v,i,j|
            vv << v
            ii << i
            jj << j
          end

          expect(vv).to eq([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 16, 17, 18, 19, 20])
          expect(ii).to eq([[0]*5, [1]*5, [2]*5, [3]*1, [4]*5].flatten)
          expect(jj).to eq([0,1,2,3,4,  0,1,2,3,5,  0,1,2,4,5,  3,  0,2,3,4,5])
        end

        it "should visit each stored element of the slice in order by indices" do

          vv = []
          ii = []
          jj = []
          m.each_ordered_stored_with_indices do |v,i,j|
            vv << v
            ii << i
            jj << j
          end
          expect(ii).to eq([0,0,0, 1,1,   2,  3,3  ])
          expect(jj).to eq([0,1,2, 0,1,   2,  1,2  ])
          expect(vv).to eq([7,8,9, 12,13, 0, 17,18 ])
        end
      end

      it "should visit each cell in the matrix as if dense, making indices available" do
        vv = []
        ii = []
        jj = []
        n.each_with_indices do |v,i,j|
          vv << v
          ii << i
          jj << j
        end

        expect(vv).to eq([1,2,3,4,5,0,6,7,8,9,0,10,11,12,13,0,14,15,0,0,0,0,0,0,16,0,17,18,19,20])
        expect(ii).to eq([[0]*6, [1]*6, [2]*6, [3]*6, [4]*6].flatten)
        expect(jj).to eq([0,1,2,3,4,5]*5)
      end

      it "should visit each cell in the slice as if dense, making indices available" do
        vv = []
        ii = []
        jj = []
        m.each_with_indices do |v,i,j|
          vv << v
          ii << i
          jj << j
        end
        expect(jj).to eq([0,1,2]*4)
        expect(ii).to eq([[0]*3, [1]*3, [2]*3, [3]*3].flatten)
        expect(vv).to eq([7,8,9,12,13,0,0,0,0,0,17,18])

      end

      if stype == :list or stype == :dense then
        it "should correctly map to a matrix with a single element" do
          nm = N.new([1], [2.0], stype: stype)
          expect(nm.map { |e| e**2 }).to eq N.new([1], [4.0], stype: stype)
        end

        it "should correctly map to a matrix with multiple elements" do
          nm = N.new([2], [2.0, 2.0], stype: stype)
          expect(nm.map { |e| e**2 }).to eq N.new([2], [4.0, 4.0], stype: stype)
        end
      end
    end
  end
end
