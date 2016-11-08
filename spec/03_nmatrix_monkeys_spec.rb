require_relative 'spec_helper'

describe NMatrix do
  describe "#to_a" do
    it "creates an Array with the same dimensions" do
      n = NMatrix.seq([3,2])
      expect(n.to_a).to eq([[0, 1], [2, 3], [4, 5]])
    end

    it "creates an Array with the proper element type" do
      n = NMatrix.seq([3,2], dtype: :float64)
      expect(n.to_a).to eq([[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]])
    end

    it "properly interprets list matrices" do
      n = NMatrix.seq([3,2], stype: :list)
      expect(n.to_a).to eq([[0, 1], [2, 3], [4, 5]])
    end

    it "properly interprets yale matrices" do
      n = NMatrix.seq([3,2], stype: :yale)
      expect(n.to_a).to eq([[0, 1], [2, 3], [4, 5]])
    end
  end
end

describe Array do
  describe "#to_nm" do
    # [0, 1, 2, 3, 4, 5]
    let(:a) {(0..5).to_a}

    it "uses a given shape and type" do
      expect(a.to_nm([3,2]).dtype).to eq :int64
      expect(a.to_nm([3,2])).to eq(NMatrix.seq([3,2]))
    end

    it "guesses dtype based on first element" do
      a[0] = 0.0
      expect(a.to_nm([3,2]).dtype).to eq :float64
    end

    it "defaults to dtype :object if necessary" do
      a = %w(this is an array of strings)
      expect(a.to_nm([3,2]).dtype).to eq :object
      expect(a.to_nm([3,2])).to eq(NMatrix.new([3,2], a, dtype: :object))
    end

    it "attempts to intuit the shape of the Array" do
      a = [[0, 1], [2, 3], [4, 5]]
      expect(a.to_nm).to eq(NMatrix.new([3,2], a.flatten))
      expect(a.to_nm.dtype).to eq :int64
    end

    it "creates an object Array for inconsistent dimensions" do
      a = [[0, 1, 2], [3], [4, 5]]
      expect(a.to_nm).to eq(NMatrix.new([3], a, dtype: :object))
      expect(a.to_nm.dtype).to eq :object
    end

    it "intuits shape of Array into multiple dimensions" do
      a = [[[0], [1]], [[2], [3]], [[4], [5]]]
      expect(a.to_nm).to eq(NMatrix.new([3,2,1], a.flatten))
      expect(a).to eq(a.to_nm.to_a)
    end

    it "is reflective with NMatrix#to_a" do
      a = [[0, 1, 2], [3], [4, 5]]
      expect(a).to eq(a.to_nm.to_a)
    end

    it "does not permanently alter the Array" do
      a = [[0, 1], [2, 3], [4, 5]]
      expect(a.to_nm).to eq(NMatrix.new([3,2], a.flatten))
      expect(a).to eq([[0, 1], [2, 3], [4, 5]])
    end
  end
end

