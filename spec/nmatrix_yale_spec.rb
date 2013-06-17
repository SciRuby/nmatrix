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
# == nmatrix_yale_spec.rb
#
# Basic tests for NMatrix's Yale storage type.
#
require "./lib/nmatrix"

describe NMatrix do
  context :yale do
    it "calculates itype" do
      NMatrix.itype_by_shape([4,4]).should == :uint8
      NMatrix.itype_by_shape(4).should == :uint8
      ## FIXME: Check larger shapes for correct itype
    end

    it "compares two empty matrices" do
      n = NMatrix.new(:yale, [4,4], :float64)
      m = NMatrix.new(:yale, [4,4], :float64)
      n.should == m
    end

    it "compares two matrices following basic assignments" do
      n = NMatrix.new(:yale, [2,2], :float64)
      m = NMatrix.new(:yale, [2,2], :float64)
      m[0,0] = 1
      m[0,1] = 1
      n.should_not == m
      n[0,0] = 1
      n.should_not == m
      n[0,1] = 1
      n.should == m
    end

    it "compares two matrices following elementwise operations" do
      n = NMatrix.new(:yale, [2,2], :float64)
      n[0,1] = 1
      m = NMatrix.new(:yale, [2,2], :float64)
      m[0,1] = -1
      r = NMatrix.new(:yale, [2,2], :float64)
      r[0,1] = 0
      (n+m).should == r
    end

    it "sets diagonal values" do
      n = NMatrix.new(:yale, [2,3], :float64)
      n.extend(NMatrix::YaleFunctions)
      n[1,1] = 0.1
      n[0,0] = 0.2
      n.yale_d.should == [0.2, 0.1]
    end

    it "gets non-diagonal rows as hashes" do
      n = NMatrix.new(:yale, [4,6], :float64)
      n.extend(NMatrix::YaleFunctions)
      n[0,0] = 0.1
      n[0,2] = 0.2
      n[0,3] = 0.3
      n[1,5] = 0.4
      h = n.yale_nd_row(0, :hash)
      h.should == {2 => 0.2, 3 => 0.3}
    end

    it "gets non-diagonal occupied column indices for a given row" do
      n = NMatrix.new(:yale, [4,6], :float64)
      n.extend(NMatrix::YaleFunctions)
      n[0,0] = 0.1
      n[0,2] = 0.2
      n[0,3] = 0.3
      n[1,5] = 0.4
      a = n.yale_nd_row(0, :array)
      a.should == [2,3]
    end

    it "does not resize until necessary" do
      n = NMatrix.new(:yale, [2,3], :float64)
      n.extend(NMatrix::YaleFunctions)
      n.yale_size.should == 3
      n.capacity.should == 5
      n[0,0] = 0.1
      n[0,1] = 0.2
      n[1,0] = 0.3
      n.yale_size.should == 5
      n.capacity.should == 5
    end


    it "sets when not resizing" do
      n = NMatrix.new(:yale, [2,3], :float64)
      n.extend(NMatrix::YaleFunctions)
      n[0,0] = 0.1
      n[0,1] = 0.2
      n[1,0] = 0.3
      n.yale_a == [0.1, 0.0, 0.0, 0.2, 0.3]
      n.yale_ija == [3,4,5,1,0]
    end

    it "sets when resizing" do
      n = NMatrix.new(:yale, [2,3], :float64)
      n.extend(NMatrix::YaleFunctions)
      n[0,0] = 0.01
      n[1,1] = 0.1
      n[0,1] = 0.2
      n[1,0] = 0.3
      n[1,2] = 0.4
      n.yale_d.should == [0.01, 0.1]
      n.yale_ia.should == [3,4,6]
      n.yale_ja.should == [1,0,2,nil]
      n.yale_lu.should == [0.2, 0.3, 0.4, nil]
    end

    it "resizes without erasing values" do
      require 'yaml'

      associations = File.open('spec/nmatrix_yale_resize_test_associations.yaml') { |y| YAML::load(y) }

      n = NMatrix.new(:yale, [618, 2801], associations.size, :byte)

      associations.each_pair do |j,i|
        n[i,j] = 1
        n[i,j].should be(1), "Value at #{i},#{j} not inserted correctly!"
      end

      associations.each_pair do |j,i|
        n[i,j].should be(1), "Value at #{i},#{j} erased during resize!"
      end
    end

    it "sets values within rows" do
      n = NMatrix.new(:yale, [3,20], :float64)
      n.extend(NMatrix::YaleFunctions)
      n[2,1]   = 1.0
      n[2,0]   = 1.5
      n[2,15]  = 2.0
      n.yale_lu.should == [1.5, 1.0, 2.0]
      n.yale_ja.should == [0, 1, 15]
    end

    it "gets values within rows" do
      n = NMatrix.new(:yale, [3,20], :float64)
      n[2,1]   = 1.0
      n[2,0]   = 1.5
      n[2,15]  = 2.0
      n[2,1].should == 1.0
      n[2,0].should == 1.5
      n[2,15].should == 2.0
    end

    it "sets values within large rows" do
      n = NMatrix.new(:yale, [10,300], :float64)
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
      n.yale_lu.should == [1.5, 1.0, 2.0, 8.0, 5.0, 6.0, 3.0, 4.0, 2.0, 7.0]
      n.yale_ja.should == [0,   1,   15,  100, 289, 290, 291, 292, 293, 299]
    end

    it "gets values within large rows" do
      n = NMatrix.new(:yale, [10,300], :float64)
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
        n[5,j].should == n.yale_lu[idx]
      end
    end

    it "dots two identical matrices" do
      a = NMatrix.new(:yale, 4, :float64)
      a[0,1] = 4.0
      a[1,2] = 1.0
      a[1,3] = 1.0
      a[3,1] = 2.0

      b = a.dup
      c = a.dot b

      c[0,0].should == 0.0
      c[0,1].should == 0.0
      c[0,2].should == 4.0
      c[0,3].should == 4.0
      c[1,0].should == 0.0
      c[1,1].should == 2.0
      c[1,2].should == 0.0
      c[1,3].should == 0.0
      c[2,0].should == 0.0
      c[2,1].should == 0.0
      c[2,2].should == 0.0
      c[2,3].should == 0.0
      c[3,0].should == 0.0
      c[3,1].should == 0.0
      c[3,2].should == 2.0
      c[3,3].should == 2.0
    end

    it "dots two identical matrices where a positive and negative partial sum cancel on the diagonal" do
      a = NMatrix.new(:yale, 4, :float64)

      a[0,0] = 1.0
      a[0,1] = 4.0
      a[1,2] = 2.0
      a[1,3] = -4.0
      a[3,1] = 4.0
      a[3,3] = 4.0

      b = a.dup
      c = a.dot b

      #c[0,0].should == 1.0
      #c[0,1].should == 4.0
      #c[0,2].should == 8.0
      #c[0,3].should == -16.0
      #c[1,0].should == 0.0
      #c[1,1].should == -16.0
      #c[1,2].should == 0.0
      #c[1,3].should == -16.0
      #c[2,0].should == 0.0
      #c[2,1].should == 0.0
      #c[2,2].should == 0.0
      #c[2,3].should == 0.0
      #c[3,0].should == 0.0
      #c[3,1].should == 0.0
      #c[3,2].should == 8.0
      #c[3,3].should == 0.0 # this is the positive and negative partial sum cancel

      c.extend(NMatrix::YaleFunctions)

      c.yale_ija.reject { |i| i.nil? }.should == [5,8,9,9,11,1,2,3,3,1,2]
      c.yale_a.reject { |i| i.nil? }.should == [1.0, -16.0, 0.0, 0.0, 0.0, 4.0, 8.0, -16.0, -16.0, 16.0, 8.0]

    end

    it "dots two vectors" do
      pending "causes segmentation fault"
      n = NVector.new(:yale, [8,1], :int64)
      m = NVector.new(:yale, [1,8], :int64)
      n.dot(m)
      m.dot(n)
    end

    it "transposes" do
      a = NMatrix.new(:yale, 4, :float64)
      a[0,0] = 1.0
      a[0,1] = 4.0
      a[1,2] = 2.0
      a[1,3] = -4.0
      a[3,1] = 5.0
      a[3,3] = 6.0
      b = a.transpose

      b[0,0].should == 1.0
      b[1,0].should == 4.0
      b[2,0].should == 0.0
      b[3,0].should == 0.0
      b[0,1].should == 0.0
      b[1,1].should == 0.0
      b[2,1].should == 2.0
      b[3,1].should == -4.0
      b[0,3].should == 0.0
      b[1,3].should == 5.0
      b[2,3].should == 0.0
      b[3,3].should == 6.0
    end
  end
end
