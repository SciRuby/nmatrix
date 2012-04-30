require "nmatrix"

describe "Slicing" do
  before :each do
    @m = NMatrix.new(:dense, [3,3], (0..9).to_a, :int32)
  end

  it 'should return NMatrix' do
    n = @m[0..1,0..1]
    @m.should == NMatrix.new(:dense, [2,2], [0,1,3,4], :int32) #== don't work
  end

  it 'should return [2x2] matrix with refs to self elements' do
    n = @m[1..2,0..1]
    n.shape.should eql([2,2])
    
    n[0,0].should == @m[1,0]
    n[0,0] = -9
    @m[1,0].should eql(-9)
  end

  it 'should return [1x2] matrix with refs to self elements' do
    n = @m[0,1..2]
    n.shape.should eql([1,2])
    
    n[0,0].should == @m[0,1]
    n[0,0] = -9
    @m[0,1].should eql(-9)
  end

  it 'should return [2x1] matrix with refs to self elements' do
    n = @m[0..1,1]
    n.shape.should eql([2,1])
    
    n[0,0].should == @m[0,1]
    n[0,0] = -9
    @m[0,1].should eql(-9)
  end

  it 'should set value from NMatrix'
    
  it 'should be elimenanted GC correctly'  do
    1.times do
      n = @m[1..2,0..1]
    end
    GC.start
    @m.should == NMatrix.new(:dense, [3,3], (0..9).to_a, :int32)

    n = nil
    1.times do
      m = NMatrix.new(:dense, [2,2], [1,2,3,4])
      n = m[0..1,0..1]
    end
    GC.start
    n.should == NMatrix.new(:dense, [2,2], [1,2,3,4])
  end
end
