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
# == stat_spec.rb
#
# Tests for statistical functions in NMatrix.
#

require_relative 'spec_helper'
require 'pry'

describe "Statistical functions" do
  context "mapping and reduction related functions" do
    [:dense, :yale, :list].each do |stype|
      context "on #{stype} matrices" do
        let(:nm_1d) { NMatrix.new([5], [5.0,0.0,1.0,2.0,3.0], stype: stype) unless stype == :yale }
        let(:nm_2d) { NMatrix.new([2,2], [0.0, 1.0, 2.0, 3.0], stype: stype) }

        it "behaves like Enumerable#reduce with no argument to reduce" do
          expect(nm_1d.reduce_along_dim(0) { |acc, el| acc + el }.to_f).to eq 11 unless stype == :yale
          expect(nm_2d.reduce_along_dim(1) { |acc, el| acc + el }).to eq NMatrix.new([2,1], [1.0, 5.0], stype: stype)
        end

        it "should calculate the mean along the specified dimension" do
          unless stype == :yale then
            puts nm_1d.mean
            expect(nm_1d.mean).to eq NMatrix.new([1], [2.2], stype: stype, dtype: :float64)
          end
          expect(nm_2d.mean).to eq NMatrix[[1.0,2.0], stype: stype]
          expect(nm_2d.mean(1)).to eq NMatrix[[0.5], [2.5], stype: stype]
        end

        it "should calculate the minimum along the specified dimension" do
          expect(nm_1d.min).to eq 0.0 unless stype == :yale
          expect(nm_2d.min).to eq NMatrix[[0.0, 1.0], stype: stype]
          expect(nm_2d.min(1)).to eq NMatrix[[0.0], [2.0], stype: stype]
        end

        it "should calculate the maximum along the specified dimension" do
          expect(nm_1d.max).to eq 5.0  unless stype == :yale
          expect(nm_2d.max).to eq NMatrix[[2.0, 3.0], stype: stype]
        end

        it "should calculate the variance along the specified dimension" do
          expect(nm_1d.variance).to eq NMatrix[3.7, stype: stype] unless stype == :yale
          expect(nm_2d.variance(1)).to eq NMatrix[[0.5], [0.5], stype: stype]
        end

        it "should calculate the sum along the specified dimension" do
          expect(nm_1d.sum).to eq NMatrix[11.0, stype: stype] unless stype == :yale
          expect(nm_2d.sum).to eq NMatrix[[2.0, 4.0], stype: stype]
        end

        it "should calculate the standard deviation along the specified dimension" do
          expect(nm_1d.std).to eq NMatrix[Math.sqrt(3.7), stype: stype] unless stype == :yale
          expect(nm_2d.std(1)).to eq NMatrix[[Math.sqrt(0.5)], [Math.sqrt(0.5)], stype: stype]
        end

        it "should raise an ArgumentError when any invalid dimension is provided" do
          expect { nm_1d.mean(3) }.to raise_exception(RangeError) unless stype == :yale
          expect { nm_2d.mean(3) }.to raise_exception(RangeError)
        end

        it "should convert to float if it contains only a single element" do
          expect(NMatrix[4.0, stype: stype].to_f).to eq 4.0  unless stype == :yale
          expect(NMatrix[[[[4.0]]], stype: stype].to_f).to eq 4.0  unless stype == :yale
          expect(NMatrix[[4.0], stype: stype].to_f).to eq 4.0
        end

        it "should raise an index error if it contains more than a single element" do
          expect { nm_1d.to_f }.to raise_error(IndexError)  unless stype == :yale
          expect { nm_2d.to_f }.to raise_error(IndexError)
        end

        it "should map a block to all elements" do
          expect(nm_1d.map { |e| e ** 2 }).to eq NMatrix[25.0,0.0,1.0,4.0,9.0, stype: stype] unless stype == :yale
          expect(nm_2d.map { |e| e ** 2 }).to eq NMatrix[[0.0,1.0],[4.0,9.0], stype: stype]
        end

        it "should map! a block to all elements in place" do
          fct = Proc.new { |e| e ** 2 }
          unless stype == :yale then
            expected1 = nm_1d.map(&fct)
            nm_1d.map!(&fct)
            expect(nm_1d).to eq expected1
          end
          expected2 = nm_2d.map(&fct)
          nm_2d.map!(&fct)
          expect(nm_2d).to eq expected2
        end

        it "should return an enumerator for map without a block" do
          expect(nm_2d.map).to be_a Enumerator
        end

        it "should return an enumerator for reduce without a block" do
          expect(nm_2d.reduce_along_dim(0)).to be_a Enumerator
        end

        it "should return an enumerator for each_along_dim without a block" do
          expect(nm_2d.each_along_dim(0)).to be_a Enumerator
        end

        it "should iterate correctly for map without a block" do
          en = nm_1d.map unless stype == :yale
          expect(en.each { |e| e**2 }).to eq nm_1d.map { |e| e**2 } unless stype == :yale
          en = nm_2d.map
          expect(en.each { |e| e**2 }).to eq nm_2d.map { |e| e**2 }
        end

        it "should iterate correctly for reduce without a block" do
          unless stype == :yale then
            en = nm_1d.reduce_along_dim(0, 1.0)
            expect(en.each { |a, e| a+e }.to_f).to eq 12
          end
          en = nm_2d.reduce_along_dim(1, 1.0)
          expect(en.each { |a, e| a+e }).to eq NMatrix[[2.0],[6.0], stype: stype]
        end

        it "should iterate correctly for each_along_dim without a block" do
          unless stype == :yale then
            res = NMatrix.zeros_like(nm_1d[0...1])
            en = nm_1d.each_along_dim(0)
            en.each { |e| res += e }
            expect(res.to_f).to eq 11
          end
          res = NMatrix.zeros_like (nm_2d[0...2, 0])
          en = nm_2d.each_along_dim(1)
          en.each { |e| res += e }
          expect(res).to eq NMatrix[[1.0], [5.0], stype: stype]
        end

        it "should yield matrices of matching dtype for each_along_dim" do
          m = NMatrix.new([2,3], [1,2,3,3,4,5], dtype: :complex128, stype: stype)
          m.each_along_dim(1) do |sub_m|
            expect(sub_m.dtype).to eq :complex128
          end
        end

        it "should reduce to a matrix of matching dtype for reduce_along_dim" do
          m = NMatrix.new([2,3], [1,2,3,3,4,5], dtype: :complex128, stype: stype)
          m.reduce_along_dim(1) do |acc, sub_m|
            expect(sub_m.dtype).to eq :complex128
            acc
          end

          m.reduce_along_dim(1, 0.0) do |acc, sub_m|
            expect(sub_m.dtype).to eq :complex128
            acc
          end
        end

        it "should allow overriding the dtype for reduce_along_dim" do
          m = NMatrix[[1,2,3], [3,4,5], dtype: :complex128]
          m.reduce_along_dim(1, 0.0, :float64) do |acc, sub_m|
            expect(acc.dtype).to eq :float64
            acc
          end

          m = NMatrix[[1,2,3], [3,4,5], dtype: :complex128, stype: stype]
          m.reduce_along_dim(1, nil, :float64) do |acc, sub_m|
            expect(acc.dtype).to eq :float64
            acc
          end
        end

        it "should convert integer dtypes to float when calculating mean" do
          m = NMatrix[[1,2,3], [3,4,5], dtype: :int32, stype: stype]
          expect(m.mean(0).dtype).to eq :float64
        end

        it "should convert integer dtypes to float when calculating variance" do
          m = NMatrix[[1,2,3], [3,4,5], dtype: :int32, stype: stype]
          expect(m.variance(0).dtype).to eq :float64
        end

        it "should convert integer dtypes to float when calculating standard deviation" do
          m = NMatrix[[1,2,3], [3,4,5], dtype: :int32, stype: stype]
          expect(m.std(0).dtype).to eq :float64
        end
      end
    end
  end
end
