require 'spec_helper'

describe 'issue 411' do
  # positive definite
  it "returns false for a non-square matrix" do 
    n = NMatrix.random([3,4])
    expect(n.positive_definite?).to eq(false)
  end

  it "returns false for a negative 1x1 matrix" do
    n = NMatrix.new([1,1], [-1])
    expect(n.positive_definite?).to eq(false)
  end

  it "returns true for a positive 1x1 matrix" do
    n = NMatrix.new([1,1], [1])
    expect(n.positive_definite?).to eq(true)
  end

  it "returns false for a square matrix with negative diagonal value" do
    n = NMatrix.new([2,2], [0, 1, 2, -3])
    expect(n.positive_definite?).to eq(false)
  end

  it "returns false for a 2x2 matrix with a negative determinant of submatrix" do
    n = NMatrix.new([2,2], [1, 10, 20, 2])
    expect(n.positive_definite?).to eq(false)
  end

  it "returns false for a 3x3 matrix with a negative determinant of submatrix" do
    n = NMatrix.new([3,3], [1,10,4,20,2,5,8,7,6])
    expect(n.positive_definite?).to eq(false)
  end

  it "returns false for a 3x3 matrix with all zero values" do
    n = NMatrix.zeros([3,3])
    expect(n.positive_definite?).to eq(false)
  end

  it "returns true for a 2x2 matrix with positive submatrix determinants" do
    n = NMatrix.new([2,2], [20, 1, 10, 2])
    expect(n.positive_definite?).to eq(true)
  end

  it "returns true for a 3x3 identity matrix" do
    n = NMatrix.eye(3)
    expect(n.positive_definite?).to eq(true)
  end

  it "returns true for a 3x3 matrix with positive submatrix determinants" do
    n = NMatrix.new([3,3], [2,-1,-1,-1,2,1,-1,1,2])
    expect(n.positive_definite?).to eq(true)
  end

  it "returns false for a 3x3 matrix with one nonpositive (0) submatrix determinant" do
    n = NMatrix.new([3,3], [2,-1,-1,-1,2,-1,-1,-1,2])
    expect(n.positive_definite?).to eq(false)
  end

  it "returns false for a 3x3 matrix with one negative submatrix determinant" do
    n = NMatrix.new([3,3], [1,2,3,2,5,4,3,4,9])
    expect(n.positive_definite?).to eq(false)
  end

  it "returns true for a 4x4 matrix with positive submatrix determinants" do
    n = NMatrix.new([4,4], [1,2,0,0, 2,6,-2, 0, 0, -2,5,-2,0,0,-2,3])
    expect(n.positive_definite?).to eq(true)
  end

  # positive semidefinite
  it "returns true for a 3x3 matrix with one nonpositive (0) submatrix determinant" do
    n = NMatrix.new([3,3], [2,-1,-1,-1,2,-1,-1,-1,2])
    expect(n.positive_semidefinite?).to eq(true)
  end

  it "returns false for a 3x3 negative definite matrix" do
    n = NMatrix.new([3,3], [-10,1,10,1,1,-10,19,20,1])
    expect(n.positive_semidefinite?).to eq(false)
  end
end
