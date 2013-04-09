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
# == slec_spec.rb
#
# Test of slice operations.
#
require File.dirname(__FILE__) + "/spec_helper.rb"

describe "Slice operation" do
  [:dense, :list, :yale].each do |stype|
    context "for #{stype}" do
      before :each do
        @m = create_matrix(stype)
      end

      unless stype == :yale
        it "should have #is_ref? method" do
          a = @m[0..1, 0..1]
          b = @m.slice(0..1, 0..1)


          @m.is_ref?.should be_false
          a.is_ref?.should be_true
          b.is_ref?.should be_false
        end

        it "reference should compare with non-reference" do
          @m.slice(1..2,0..1).should == @m[1..2, 0..1]
          @m[1..2,0..1].should == @m.slice(1..2, 0..1)
          @m[1..2,0..1].should == @m[1..2, 0..1]
        end
      end # unless stype == :yale

      context "with copying" do
        it 'should return an NMatrix' do
          n = @m.slice(0..1,0..1)
          nm_eql(n, NMatrix.new([2,2], [0,1,3,4], :int32)).should be_true
        end

        it 'should return a copy of 2x2 matrix to self elements' do
          n = @m.slice(1..2,0..1)
          n.shape.should eql([2,2])

          n[1,1].should == @m[2,1]
          n[1,1] = -9
          @m[2,1].should eql(7)
        end

        it 'should return a 1x2 matrix with refs to self elements' do
          n = @m.slice(0,1..2)
          n.shape.should eql([1,2])

          n[0,0].should == @m[0,1]
          n[0,0] = -9
          @m[0,1].should eql(1)
        end

        it 'should return a 2x1 matrix with refs to self elements' do
          n = @m.slice(0..1,1)
          n.shape.should eql([2,1])

          n[0,0].should == @m[0,1]
          n[0,0] = -9
          @m[0,1].should eql(1)
        end

        it 'should be correct slice for range 0..2 and 0...3' do
          @m.slice(0..2,0..2).should == @m.slice(0...3,0...3)
        end

        [:dense, :list, :yale].each do |cast_type|
          it "should cast from #{stype.upcase} to #{cast_type.upcase}" do
            nm_eql(@m.slice(1..2, 1..2).cast(cast_type, :int32), @m.slice(1..2,1..2)).should be_true
            nm_eql(@m.slice(0..1, 1..2).cast(cast_type, :int32), @m.slice(0..1,1..2)).should be_true
            nm_eql(@m.slice(1..2, 0..1).cast(cast_type, :int32), @m.slice(1..2,0..1)).should be_true
            nm_eql(@m.slice(0..1, 0..1).cast(cast_type, :int32), @m.slice(0..1,0..1)).should be_true

            # Non square
            nm_eql(@m.slice(0..2, 1..2).cast(cast_type, :int32), @m.slice(0..2,1..2)).should be_true
            nm_eql(@m.slice(1..2, 0..2).cast(cast_type, :int32), @m.slice(1..2,0..2)).should be_true

            # Full
            nm_eql(@m.slice(0..2, 0..2).cast(cast_type, :int32), @m).should be_true
          end
        end
      end

      if stype == :yale
        context "by reference" do
          it "should raise an error" do
            expect{ @m[1..2,1..2] }.to raise_error(NotImplementedError)
          end
        end

        context "by copy" do
          it "should correctly preserve zeros" do
            @m = NMatrix.new(:yale, 3, :int64)
            column_slice = @m.column(2, :copy)
            column_slice[0,0].should == 0
            column_slice[1,0].should == 0
            column_slice[2,0].should == 0
          end
        end
      else
        context "by reference" do
          it 'should return an NMatrix' do
            n = @m[0..1,0..1]
            nm_eql(n, NMatrix.new([2,2], [0,1,3,4], :int32)).should be_true
          end

          it 'should return a 2x2 matrix with refs to self elements' do
            n = @m[1..2,0..1]
            n.shape.should eql([2,2])

            n[0,0].should == @m[1,0]
            n[0,0] = -9
            @m[1,0].should eql(-9)
          end

          it 'should return a 1x2 matrix with refs to self elements' do
            n = @m[0,1..2]
            n.shape.should eql([1,2])

            n[0,0].should == @m[0,1]
            n[0,0] = -9
            @m[0,1].should eql(-9)
          end

          it 'should return a 2x1 matrix with refs to self elements' do
            n = @m[0..1,1]
            n.shape.should eql([2,1])

            n[0,0].should == @m[0,1]
            n[0,0] = -9
            @m[0,1].should eql(-9)
          end

          it 'should set value from NMatrix'

          it 'should slice again' do
            n = @m[1..2, 1..2]
            nm_eql(n[1,0..1], NMatrix.new([1,2], [7,8])).should be_true
          end

          it 'should be correct slice for range 0..2 and 0...3' do
            @m[0..2,0..2].should == @m[0...3,0...3]
          end

          if stype == :dense
            [:byte,:int8,:int16,:int32,:int64,:float32,:float64,:rational64,:rational128].each do |left_dtype|
              [:byte,:int8,:int16,:int32,:int64,:float32,:float64,:rational64,:rational128].each do |right_dtype|

                # Won't work if they're both 1-byte, due to overflow.
                next if [:byte,:int8].include?(left_dtype) && [:byte,:int8].include?(right_dtype)

                # For now, don't bother testing int-int mult.
                #next if [:int8,:int16,:int32,:int64].include?(left_dtype) && [:int8,:int16,:int32,:int64].include?(right_dtype)
                it "handles #{left_dtype.to_s} dot #{right_dtype.to_s} matrix multiplication" do
                  #STDERR.puts "dtype=#{dtype.to_s}"
                  #STDERR.puts "2"

                  nary = if left_dtype.to_s =~ /complex/
                           COMPLEX_MATRIX43A_ARRAY
                         elsif left_dtype.to_s =~ /rational/
                           RATIONAL_MATRIX43A_ARRAY
                         else
                           MATRIX43A_ARRAY
                         end

                  mary = if right_dtype.to_s =~ /complex/
                           COMPLEX_MATRIX32A_ARRAY
                         elsif right_dtype.to_s =~ /rational/
                           RATIONAL_MATRIX32A_ARRAY
                         else
                           MATRIX32A_ARRAY
                         end

                  n = NMatrix.new([4,3], nary, left_dtype)[1..3,1..2]
                  m = NMatrix.new([3,2], mary, right_dtype)[1..2,0..1]

                  r = n.dot m
                  r.shape.should eql([3,2])

                  r[0,0].should == 219.0
                  r[0,1].should == 185.0
                  r[1,0].should == 244.0
                  r[1,1].should == 205.0
                  r[2,0].should == 42.0
                  r[2,1].should == 35.0

                end
              end
            end
          end

          it 'should be cleaned up by garbage collector without errors'  do
            1.times do
              n = @m[1..2,0..1]
            end
            GC.start
            @m.should == NMatrix.new(:dense, [3,3], (0..9).to_a, :int32).cast(stype, :int32)
            n = nil
            1.times do
              m = NMatrix.new(:dense, [2,2], [1,2,3,4]).cast(stype, :int32)
              n = m[0..1,0..1]
            end
            GC.start
            n.should == NMatrix.new(:dense, [2,2], [1,2,3,4]).cast(stype, :int32)
          end

          [:dense, :list, :yale].each do |cast_type|
            it "should cast from #{stype.upcase} to #{cast_type.upcase}" do
              nm_eql(@m[1..2, 1..2].cast(cast_type, :int32), @m[1..2,1..2]).should be_true
              nm_eql(@m[0..1, 1..2].cast(cast_type, :int32), @m[0..1,1..2]).should be_true
              nm_eql(@m[1..2, 0..1].cast(cast_type, :int32), @m[1..2,0..1]).should be_true
              nm_eql(@m[0..1, 0..1].cast(cast_type, :int32), @m[0..1,0..1]).should be_true

              # Non square
              nm_eql(@m[0..2, 1..2].cast(cast_type, :int32), @m[0..2,1..2]).should be_true
              nm_eql(@m[1..2, 0..2].cast(cast_type, :int32), @m[1..2,0..2]).should be_true

              # Full
              nm_eql(@m[0..2, 0..2].cast(cast_type, :int32), @m).should be_true
            end
          end
        end
      end
    end # unless stype == :yale
  end

  # Stupid but independent comparison
  def nm_eql(n, m) #:nodoc:
    if n.shape != m.shape
      false
    else
      n.shape[0].times do |i|
        n.shape[1].times do |j|
          if n[i,j] != m[i,j]
            puts "n[#{i},#{j}] != m[#{i},#{j}] (#{n[i,j]} != #{m[i,j]})"
            return false
          end
        end
      end
      true
    end
  end
end
