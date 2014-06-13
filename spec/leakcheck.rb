require "./lib/nmatrix"

# Fixed:
#n = NMatrix.new(:yale, [8,2], :int64)
#m = NMatrix.new(:yale, [2,8], :int64)
#100.times do
#  n.dot(m)
#end
#GC.start

# Remaining:
100.times do |t|
  n = NMatrix.new(:dense, 1000, :float64)
  n[0,t] = 1.0
  puts n[t,0]
end
