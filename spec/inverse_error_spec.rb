require File.join(File.dirname(__FILE__), "spec_helper.rb")
describe 'This causes a memory corruption for me' do 
  a = NMatrix.new([5,6], 
    %w|8.79 9.93 9.83 5.45 3.16
      6.11 6.91 5.04 -0.27 7.98 
      -9.15 -7.93 4.86 4.85 3.01 
      9.57 1.64 8.83 0.74 5.80 
      -3.49 4.02 9.80 10.00 4.27 
      9.84 0.15 -8.99 -6.02 -5.31|.map(&:to_f), 
  :float32)
  a.transpose.inverse
end
