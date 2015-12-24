require_relative 'lib/nmatrix/lapacke'
require_relative 'lib/nmatrix/math'
 
# this works: n = NMatrix.new([4,3], [1,2,3, 81,51,89, 30, 2, 86, 92, 89, 60], :dtype => :float64)
n = NMatrix.new([3,3], [1,2,3,3,6,9,92,89,60], :dtype => :float32)

Float::EPSILON = 0.00000011920928955078125

sigmas = n.gesvd[1].to_a.flatten
p "sigmas: #{sigmas}"

tol = n.shape.max * sigmas.max * Float::EPSILON 
p "tolerance: #{tol}"

rank = sigmas.map { |x| x > tol ? 1 : 0 }.reduce(:+)
p "rank: #{rank}"
