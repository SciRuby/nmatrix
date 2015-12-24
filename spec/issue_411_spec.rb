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

describe "rank of matrices.  requires lapacke" do
  [:float32, :float64, :complex64, :complex128].each do |dtype|
    context dtype do
      # full rank
      it "returns full rank for identity matrix" do
        n = NMatrix.new([3,3], [1,0,0,0,1,0,0,0,1], :dtype => dtype)

        begin 
          expect(n.full_rank?).to be(true)
        rescue NotImplementedError => e
          pending e.to_s
        end
      end

      it "returns full rank for a full rank square matrix" do
        n = NMatrix.new([3,3], [81,51,89, 30, 2, 86, 92, 89, 60], :dtype => dtype)

        begin 
          expect(n.full_rank?).to be(true)
        rescue NotImplementedError => e
          pending e.to_s
        end
      end

      it "returns full rank for a full rank rectangular matrix" do
        n = NMatrix.new([4,3], [1,2,3, 81,51,89, 30, 2, 86, 92, 89, 60], :dtype => dtype)

        begin 
          expect(n.full_rank?).to be(true)
        rescue NotImplementedError => e
          pending e.to_s
        end
      end

      # rank deficient
      it "returns rank deficient for a non-full rank matrix" do
        n = NMatrix.new([3,3], [1,2,3,3,6,9,6,12,18], :dtype => dtype)

        begin 
          expect(n.rank_deficient?).to be(true)
        rescue NotImplementedError => e
          pending e.to_s
        end
      end

      # get rank of a matrix
      it "returns rank 3 for a matrix of rank 3" do
        n = NMatrix.new([3,3], [81,51,89, 30, 2, 86, 92, 89, 60], :dtype => dtype)

        begin 
          expect(n.get_rank).to be(3)
        rescue NotImplementedError => e
          pending e.to_s
        end
      end

      it "returns rank 2 for a matrix with 2 independent rows" do
        n = NMatrix.new([3,3], [1,2,3,3,6,9, 92, 89, 60], :dtype => dtype)

        begin 
          expect(n.get_rank).to be(2)
        rescue NotImplementedError => e
          pending e.to_s
        end
      end

      it "returns rank 1 for a matrix with a matrix with linearly dependent rows" do
        n = NMatrix.new([3,3], [1,2,3,3,6,9,6,12,18], :dtype => dtype)

        begin 
          expect(n.get_rank).to be(1)
        rescue NotImplementedError => e
          pending e.to_s
        end
      end
    end
  end
end
