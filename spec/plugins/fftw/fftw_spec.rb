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
# == fftw_spec.rb
#
# Tests for interfaces that are only exposed by nmatrix-fftw
#

require 'spec_helper'
require "./lib/nmatrix/fftw"

describe NMatrix do
  context "#fft" do
    before do
      @answer = NMatrix.new([10],
        [ 
          Complex(330.3200,0.0000)   , Complex(-8.4039  ,-150.3269),
          Complex(-99.4807,-68.6579) , Complex(-143.6861, -20.4273),
          Complex(67.6207  ,  8.5236), Complex(130.7800 ,  0.0000),
          Complex(67.6207 ,  -8.5236), Complex(-143.6861, 20.4273),
          Complex(-99.4807 , 68.6579), Complex(-8.4039  ,150.3269)
        ], dtype: :complex128)      
    end

    it "computes an FFT of a complex NMatrix" do
      nm = NMatrix.new([10],
        [
          Complex(9.32,0), Complex(44,0), Complex(125,0), Complex(34,0),
          Complex(31,0),   Complex(44,0), Complex(12,0),  Complex(1,0),
          Complex(53.23,0),Complex(-23.23,0)], dtype: :complex128)
      expect(nm.fft.round(4)).to eq(@answer)
    end
  end

  context "#fft2" do
    it "computes 2D FFT if NMatrix has such shape" do
      input = NMatrix.new([2,2],
        [
          Complex(9.3200,0), Complex(43.0000,0),
          Complex(3.2000,0), Complex(4.0000,0)
        ], dtype: :complex128
      )
      output = NMatrix.new([2,2],
        [
          Complex(59.520,0), Complex(-34.480,0),
          Complex(45.120,0),  Complex(-32.880,0),
        ], dtype: :complex128
      )
      expect(input.fft2.round(4)).to eq(output)   
    end
  end
end

describe NMatrix::FFTW, focus: true do
  describe NMatrix::FFTW::Plan do
    context ".new" do
      it "creates a new plan for default DFT (complex input/complex output)" do
        plan = NMatrix::FFTW::Plan.new(4)
        # TODO: Figure a way to test internal C data structures.

        expect(plan.shape)    .to eq([4])
        expect(plan.size)     .to eq(4)
        expect(plan.dim)      .to eq(1)
        expect(plan.flags)     .to eq([:estimate])
        expect(plan.direction).to eq(:forward)
      end

      it "creates a new plan for multi dimensional DFT with options" do
        plan = NMatrix::FFTW::Plan.new([10,5,8],
          direction: :backward, flags: [:exhaustive, :estimate], dim: 3)

        expect(plan.shape)    .to eq([10,5,8])
        expect(plan.size)     .to eq(10*5*8)
        expect(plan.dim)      .to eq(3)
        expect(plan.flags)    .to eq([:exhaustive, :estimate])
        expect(plan.direction).to eq(:backward)
      end

      it "creates a new plan for real input/complex output" do
        plan = NMatrix::FFTW::Plan.new([5,20,10,4,2],
          direction: :forward, flags: [:patient, :exhaustive], dim: 5, 
          type: :real_complex)

        expect(plan.shape) .to eq([5,20,10,4,2])
        expect(plan.size)  .to eq(5*20*10*4*2)
        expect(plan.dim)   .to eq(5)
        expect(plan.flags) .to eq([:patient, :exhaustive])
        expect(plan.type)  .to eq(:real_complex)
      end

      it "raises error for plan with incompatible shape and dimension" do
        expect {
          NMatrix::FFTW::Plan.new([9], dim: 2, type: :real_complex)
        }.to raise_error(ArgumentError)
      end

      it "creates a new plan for real input/real output" do
        plan = NMatrix::FFTW::Plan.new([30,30], type: :real_real, 
          real_real_kind: [:rodft00, :redft10], dim: 2)

        expect(plan.shape).to eq([30,30])
        expect(plan.size) .to eq(30*30)
        expect(plan.dim)  .to eq(2)
        expect(plan.flags).to eq([:estimate])
        expect(plan.type) .to eq(:real_real)
      end

      it "creates a new plan for complex input/real output" do
        plan = NMatrix::FFTW::Plan.new([30,400], type: :complex_real, 
          dim: 2, flags: [:patient, :exhaustive])

        expect(plan.shape).to eq([30,400])
        expect(plan.size) .to eq(30*400)
        expect(plan.dim)  .to eq(2)
        expect(plan.flags).to eq([:patient, :exhaustive])
        expect(plan.type) .to eq(:complex_real)
      end
    end

    context "#set_input" do
      it "accepts nothing but complex128 input for the default or complex_real plan" do
        plan  = NMatrix::FFTW::Plan.new(4)
        input = NMatrix.new([4], [23.54,52.34,52.345,64], dtype: :float64)
        expect {
          plan.set_input(input)
        }.to raise_error(ArgumentError)

        plan = NMatrix::FFTW::Plan.new(4, type: :complex_real)
        expect {
          plan.set_input input
        }.to raise_error(ArgumentError)
      end

      it "accepts nothing but float64 input for real_complex or real_real plan" do
        plan = NMatrix::FFTW::Plan.new(4, type: :real_complex)
        input = NMatrix.new([4], [1,2,3,4], dtype: :int32)

        expect {
          plan.set_input(input)
        }.to raise_error(ArgumentError)
      end
    end

    context "#execute" do
      it "calculates a basic 1D DFT" do
        input = NMatrix.new([10],
          [
            Complex(9.32,0),
            Complex(44,0),
            Complex(125,0),
            Complex(34,0),
            Complex(31,0),
            Complex(44,0),
            Complex(12,0),
            Complex(1,0),
            Complex(53.23,0),
            Complex(-23.23,0),
          ], dtype: :complex128)

        output = NMatrix.new([10],
          [
            Complex(330.3200,0.0000),
            Complex(-8.4039  ,-150.3269),
            Complex(-99.4807 , -68.6579),
            Complex(-143.6861, -20.4273),
            Complex(67.6207  ,  8.5236),
            Complex(130.7800 ,  0.0000),
            Complex(67.6207  ,  -8.5236),
            Complex(-143.6861, 20.4273),
            Complex(-99.4807 , 68.6579),
            Complex(-8.4039  ,150.3269)
          ], dtype: :complex128)

        plan = NMatrix::FFTW::Plan.new(10)
        plan.set_input input
        expect(plan.execute).to eq(true)
        expect(plan.output.round(4)).to eq(output)
      end

      it "calculates 2D DFT with options" do
        input = NMatrix.new([2,2],
          [
            Complex(9.3200,0), Complex(43.0000,0),
            Complex(3.2000,0), Complex(4.0000,0)
          ], dtype: :complex128
        )

        output = NMatrix.new([2,2],
          [
            Complex(59.520,0), Complex(-34.480,0),
            Complex(45.120,0),  Complex(-32.880,0),
          ], dtype: :complex128
        )

        plan = NMatrix::FFTW::Plan.new([2,2],
          direction: :forward, flags: :estimate, dim: 2)
        plan.set_input input
        expect(plan.execute).to eq(true)
        expect(plan.output).to eq(output)
      end

      it "calculates ND DFT with options" do

      end

      it "calculates 1D real input/complex output DFT" do
        input  = NMatrix.new([4], [3.10, 1.73, 1.04, 2.83], dtype: :float64)
        output = NMatrix.new([3], 
          [Complex(8.70, 0), Complex(2.06, 1.1), Complex(-0.42, 0)], dtype: :complex128)
        plan = NMatrix::FFTW::Plan.new([4], type: :real_complex)
        plan.set_input input
        expect(plan.execute).to eq(true)
        expect(plan.output).to eq(output)
      end

      it "calculates 2D real input/complex output DFT" do
        input = NMatrix.new([16], [
          1  ,   5,54    ,656,
          4.3,1.32,-43.34,14 ,
          1  ,   5,    54,656,
          4.3,1.32,-43.34,14
          ], dtype: :float64) 
        output = NMatrix.new([9],
          [
            Complex(1384.56, 0.0),
            Complex(-10.719999999999999, 1327.36),
            Complex(-1320.72, 0.0),
            Complex(0.0, 0.0),
            Complex(0.0, 0.0),
            Complex(0.0, 0.0),
            Complex(1479.44, 0.0),
            Complex(-201.28, 1276.64),
            Complex(-1103.28, 0.0)
          ], dtype: :complex128
        )

        plan = NMatrix::FFTW::Plan.new([4,4], type: :real_complex, dim: 2)
        plan.set_input input
        expect(plan.execute).to eq(true)
        expect(plan.output).to eq(output)
      end

      it "calculates 1D complex input/real output DFT" do
        input = NMatrix.new([8],
          [
            Complex(9.32,0),
            Complex(43.0,0),
            Complex(3.20,0),
            Complex(4.00,0),
            Complex(5.32,0),
            Complex(3.20,0),
            Complex(4.00,0),
            Complex(5.32,0)
          ], dtype: :complex128)

        output = NMatrix.new([8], [
            115.04,59.1543,8.24,-51.1543,-72.96,-51.1543,8.24,59.1543
          ], dtype: :float64)

        plan = NMatrix::FFTW::Plan.new([8], type: :complex_real)
        plan.set_input input
        expect(plan.execute).to eq(true)
        expect(plan.output.round(4)).to eq(output)
      end

      it "calculates 2D complex input/real output DFT" do
        input = NMatrix.new([9],
          [
            Complex(9.32,0),
            Complex(43.0,0),
            Complex(3.20,0),
            Complex(4.00,0),
            Complex(5.32,0),
            Complex(3.20,0),
            Complex(4.00,0),
            Complex(5.32,0),
            Complex(45.32,0)
          ], dtype: :complex128)
        output = NMatrix.new([9], [
            118.24,-32.36,-32.36,83.86,-35.54,-33.14,83.86,-33.14,-35.54
          ], dtype: :float64)

        plan = NMatrix::FFTW::Plan.new([3,3], type: :complex_real, dim: 2)
        plan.set_input input
        expect(plan.execute).to eq(true)
        expect(plan.output.round(2)) .to eq(output)
      end

      it "calculates basic 1D real input/real output DFT of kind RODFT00" do
        input = NMatrix.new([9],
          [9.32,43.00,3.20,4.00,5.32,3.20,4.00,5.32,45.32], dtype: :float64)
        output = NMatrix.new([9],
          [126.56,28.77,165.67,-24.76,105.52,-110.31,-1.23,-116.45,-14.44],
          dtype: :float64)
        plan = NMatrix::FFTW::Plan.new([9], type: :real_real, real_real_kind: [:rodft00])
        plan.set_input input
        expect(plan.execute).to eq(true)
        expect(plan.output.round(2)).to eq(output)
      end

      it "calculates basic 1D real input/real output DFT of kind REDFT10" do
        input = NMatrix.new([9],
          [9.32,43.00,3.20,4.00,5.32,3.20,4.00,5.32,45.32], dtype: :float64)
        output = NMatrix.new([9],
          [245.36,-6.12,126.84,-62.35,35.00,-109.42,-38.24,-92.49,-21.20], 
          dtype: :float64)

        plan = NMatrix::FFTW::Plan.new([9], type: :real_real, real_real_kind: [:redft10])
        plan.set_input input
        expect(plan.execute).to eq(true)
        expect(plan.output.round(2)).to eq(output)
      end

      it "calculates 2D DFT for real input/real output of kind REDFT10, REDFT11" do
        input = NMatrix.new([9],
          [9.32,43.00,3.20,4.00,5.32,3.20,4.00,5.32,45.32], dtype: :float64)
        output = NMatrix.new([9],
          [272.181,-249.015,66.045,72.334,23.907,-228.463,85.368,-105.331,30.836],
          dtype: :float64)

        plan = NMatrix::FFTW::Plan.new([3,3], type: :real_real, 
          real_real_kind: [:redft10, :redft11], dim: 2)
        plan.set_input input
        expect(plan.execute).to eq(true)
        expect(plan.output.round(3)) .to eq(output)
      end
    end
  end
end
