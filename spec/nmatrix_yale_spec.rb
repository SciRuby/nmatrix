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
# == nmatrix_yale_spec.rb
#
# Basic tests for NMatrix's Yale storage type.
#
require 'spec_helper'
require "./lib/nmatrix"

describe NMatrix do
  context :yale do

    it "compares two empty matrices" do
      n = NMatrix.new(4, stype: :yale, dtype: :float64)
      m = NMatrix.new(4, stype: :yale, dtype: :float64)
      expect(n).to eq(m)
    end

    it "compares two matrices following basic assignments" do
      n = NMatrix.new(2, stype: :yale, dtype: :float64)
      m = NMatrix.new(2, stype: :yale, dtype: :float64)

      m[0,0] = 1
      m[0,1] = 1
      expect(n).not_to eq(m)
      n[0,0] = 1
      expect(n).not_to eq(m)
      n[0,1] = 1
      expect(n).to eq(m)
    end

    it "compares two matrices following elementwise operations" do
      n = NMatrix.new(2, stype: :yale, dtype: :float64)
      m = NMatrix.new(2, stype: :yale, dtype: :float64)
      n[0,1] = 1
      m[0,1] = -1
      x = n+m
      expect(n+m).to eq(NMatrix.new(2, 0.0, stype: :yale))
    end

    it "sets diagonal values" do
      n = NMatrix.new([2,3], stype: :yale, dtype: :float64)
      n.extend(NMatrix::YaleFunctions)
      n[1,1] = 0.1
      n[0,0] = 0.2
      expect(n.yale_d).to eq([0.2, 0.1])
    end

    it "gets non-diagonal rows as hashes" do
      n = NMatrix.new([4,6], stype: :yale, dtype: :float64)
      n.extend(NMatrix::YaleFunctions)
      n[0,0] = 0.1
      n[0,2] = 0.2
      n[0,3] = 0.3
      n[1,5] = 0.4
      h = n.yale_nd_row(0, :hash)
      expect(h).to eq({2 => 0.2, 3 => 0.3})
    end

    it "gets non-diagonal occupied column indices for a given row" do
      n = NMatrix.new([4,6], stype: :yale, dtype: :float64)
      n.extend(NMatrix::YaleFunctions)
      n[0,0] = 0.1
      n[0,2] = 0.2
      n[0,3] = 0.3
      n[1,5] = 0.4
      a = n.yale_nd_row(0, :array)
      expect(a).to eq([2,3])
    end

    it "does not resize until necessary" do
      n = NMatrix.new([2,3], stype: :yale, dtype: :float64)
      n.extend(NMatrix::YaleFunctions)
      expect(n.yale_size).to eq(3)
      expect(n.capacity).to eq(5)
      n[0,0] = 0.1
      n[0,1] = 0.2
      n[1,0] = 0.3
      expect(n.yale_size).to eq(5)
      expect(n.capacity).to eq(5)
    end


    it "sets when not resizing" do
      n = NMatrix.new([2,3], stype: :yale, dtype: :float64)
      n.extend(NMatrix::YaleFunctions)
      n[0,0] = 0.1
      n[0,1] = 0.2
      n[1,0] = 0.3
      expect(n.yale_a).to eq([0.1, 0.0, 0.0, 0.2, 0.3])
      expect(n.yale_ija).to eq([3,4,5,1,0])
    end

    it "sets when resizing" do
      n = NMatrix.new([2,3], stype: :yale, dtype: :float64)
      n.extend(NMatrix::YaleFunctions)
      n[0,0] = 0.01
      n[1,1] = 0.1
      n[0,1] = 0.2
      n[1,0] = 0.3
      n[1,2] = 0.4
      expect(n.yale_d).to eq([0.01, 0.1])
      expect(n.yale_ia).to eq([3,4,6])
      expect(n.yale_ja).to eq([1,0,2,nil])
      expect(n.yale_lu).to eq([0.2, 0.3, 0.4, nil])
    end

    it "resizes without erasing values" do
      require 'yaml'

      associations = File.open('spec/nmatrix_yale_resize_test_associations.yaml') { |y| YAML::load(y) }

      n = NMatrix.new([618,2801], stype: :yale, dtype: :byte, capacity: associations.size)
      #n = NMatrix.new(:yale, [618, 2801], associations.size, :byte)

      associations.each_pair do |j,i|
        n[i,j] = 1
        expect(n[i,j]).to be(1), "Value at #{i},#{j} not inserted correctly!"
      end

      associations.each_pair do |j,i|
        expect(n[i,j]).to be(1), "Value at #{i},#{j} erased during resize!"
      end
    end

    it "sets values within rows" do
      n = NMatrix.new([3,20], stype: :yale, dtype: :float64)
      n.extend(NMatrix::YaleFunctions)
      n[2,1]   = 1.0
      n[2,0]   = 1.5
      n[2,15]  = 2.0
      expect(n.yale_lu).to eq([1.5, 1.0, 2.0])
      expect(n.yale_ja).to eq([0, 1, 15])
    end

    it "gets values within rows" do
      n = NMatrix.new([3,20], stype: :yale, dtype: :float64)
      n[2,1]   = 1.0
      n[2,0]   = 1.5
      n[2,15]  = 2.0
      expect(n[2,1]).to eq(1.0)
      expect(n[2,0]).to eq(1.5)
      expect(n[2,15]).to eq(2.0)
    end

    it "sets values within large rows" do
      n = NMatrix.new([10,300], stype: :yale, dtype: :float64)
      n.extend(NMatrix::YaleFunctions)
      n[5,1]   = 1.0
      n[5,0]   = 1.5
      n[5,15]  = 2.0
      n[5,291] = 3.0
      n[5,292] = 4.0
      n[5,289] = 5.0
      n[5,290] = 6.0
      n[5,293] = 2.0
      n[5,299] = 7.0
      n[5,100] = 8.0
      expect(n.yale_lu).to eq([1.5, 1.0, 2.0, 8.0, 5.0, 6.0, 3.0, 4.0, 2.0, 7.0])
      expect(n.yale_ja).to eq([0,   1,   15,  100, 289, 290, 291, 292, 293, 299])
    end

    it "gets values within large rows" do
      n = NMatrix.new([10,300], stype: :yale, dtype: :float64)
      n.extend(NMatrix::YaleFunctions)
      n[5,1]   = 1.0
      n[5,0]   = 1.5
      n[5,15]  = 2.0
      n[5,291] = 3.0
      n[5,292] = 4.0
      n[5,289] = 5.0
      n[5,290] = 6.0
      n[5,293] = 2.0
      n[5,299] = 7.0
      n[5,100] = 8.0

      n.yale_ja.each_index do |idx|
        j = n.yale_ja[idx]
        expect(n[5,j]).to eq(n.yale_lu[idx])
      end
    end

    it "dots two identical matrices" do
      a = NMatrix.new(4, stype: :yale, dtype: :float64)
      a[0,1] = 4.0
      a[1,2] = 1.0
      a[1,3] = 1.0
      a[3,1] = 2.0

      b = a.dup
      c = a.dot b

      d = NMatrix.new(4, [0,0,4,4, 0,2,0,0, 0,0,0,0, 0,0,2,2], dtype: :float64, stype: :yale)

      expect(c).to eq(d)
    end

    it "dots two identical matrices where a positive and negative partial sum cancel on the diagonal" do
      a = NMatrix.new(4, 0.0, stype: :yale)

      a[0,0] = 1.0
      a[0,1] = 4.0
      a[1,2] = 2.0
      a[1,3] = -4.0
      a[3,1] = 4.0
      a[3,3] = 4.0

      b = a.dup
      c = a.dot b

      c.extend(NMatrix::YaleFunctions)

      expect(c.yale_ija.reject { |i| i.nil? }).to eq([5,8,9,9,11,1,2,3,3,1,2])
      expect(c.yale_a.reject { |i| i.nil? }).to eq([1.0, -16.0, 0.0, 0.0, 0.0, 4.0, 8.0, -16.0, -16.0, 16.0, 8.0])

    end

    it "dots two vectors" do
      n = NMatrix.new([16,1], 0, stype: :yale)
      m = NMatrix.new([1,16], 0, stype: :yale)

      n[0] = m[0] = 1
      n[1] = m[1] = 2
      n[2] = m[2] = 3
      n[3] = m[3] = 4
      n[4] = m[4] = 5
      n[5] = m[5] = 6
      n[6] = m[6] = 7
      n[7] = m[7] = 8
      n[8] = m[8] = 9
      n[15] = m[15] = 16

      nm = n.dot(m)

      # Perform the same multiplication with dense
      nmr = n.cast(:dense, :int64).dot(m.cast(:dense, :int64)).cast(:yale, :int64)

      nm.extend(NMatrix::YaleFunctions)
      nmr.extend(NMatrix::YaleFunctions)

      # We want to do a structure comparison to ensure multiplication is occurring properly, but more importantly, to
      # ensure that insertion sort is occurring as it should. If the row has more than four entries, it'll run quicksort
      # instead. Quicksort calls insertion sort for small rows, so we test both with this particular multiplication.
      expect(nm.yale_ija[0...107]).to eq(nmr.yale_ija[0...107])
      expect(nm.yale_a[0...107]).to   eq(nmr.yale_a[0...107])

      mn = m.dot(n)
      expect(mn[0,0]).to eq(541)
    end

    it "calculates the row key intersections of two matrices" do
      a = NMatrix.new([3,9], [0,1], stype: :yale, dtype: :byte, default: 0)
      b = NMatrix.new([3,9], [0,0,1,0,1], stype: :yale, dtype: :byte, default: 0)
      a.extend NMatrix::YaleFunctions
      b.extend NMatrix::YaleFunctions

      (0...3).each do |ai|
        (0...3).each do |bi|
          STDERR.puts (a.yale_ja_d_keys_at(ai) & b.yale_ja_d_keys_at(bi)).inspect
          expect(a.yale_ja_d_keys_at(ai) & b.yale_ja_d_keys_at(bi)).to eq(a.yale_row_keys_intersection(ai, b, bi))
        end
      end

    end
  end
end
