class NMatrix

  #######################################
  #
  # Shortcuts to initialize common matrices.
  #
  # by Daniel Carrera <dcarrera@hush.com>
  #
  ##########################################

  # zeros() or zeroes()
  #    
  #     Creates a new matrix of zeros with the dimensions supplied as
  #     parameters. Optional parameters include:
  #     
  #     * A storage type as the first parameter (default is :dense).
  #     * A dtype as the last parameter (default is :float64).
  #     
  # Examples:
  #   
  #   zeros(3) # =>  0.0   0.0   0.0
  #
  #   zeros(2,3,:int32) # =>  0  0  0
  #                           0  0  0
  #   
  #   zeros(:list, 5, :int32) # =>  0  0  0  0  0
  #

  def self.zeros(*params)
  	dtype = params.last.is_a?(Symbol) ? params.pop   : :float64
  	store = params.first.is_a?(Symbol) ? params.shift : :dense
	
  	NMatrix.new(store, params, 0, dtype)
  end

  alias :zeroes :zeros

  # ones()
  #    
  #     Creates a :dense matrix of ones with the dimensions supplied
  #     as parameters. Optionaly, one can specify a dtype as the last
  #     parameter (default is :float64).
  #     
  # Examples:
  #   
  #   ones(3) # =>  1.0   1.0   1.0
  #
  #   ones(2,3,:int32) # =>  1  1  1
  #                          1  1  1
  #

  def self.ones(*params)
      dtype = params.last.is_a?(Symbol) ? params.pop   : :float64
    
      NMatrix.new(params, 1, dtype)
  end

  # identity() or eye()
  #      
  #     Creates an identitiy matrix (square matrix rank 2) of the size
  #     supplied as a parameter. Optional parameters include:
  #      
  #     * A storage type as the first parameter (default is :dense).
  #     * A dtype as the last parameter (default is :float64).
  #
  # Examples:
  #      
  #    eye(3) # =>   1.0   0.0   0.0
  #                  0.0   1.0   0.0
  #                  0.0   0.0   1.0
  #    
  #    eye(3, :int32) # =>   1   0   0
  #                          0   1   0
  #                          0   0   1
  #    
  #    eye(:yale, 2, :int32) # =>   1   0
  #                                 0   1
  #

  def self.eye(*params)
  	dtype = params.last.is_a?(Symbol) ? params.pop   : :float64
  	store = params.first.is_a?(Symbol) ? params.shift : :dense
	
  	n = params[0]
  	m = zeros(store,n,n,dtype)
  	(0..n-1).each { |i| m[i,i] = 1 }
  	m
  end

  alias :identity :eye

  # seq()
  #      
  #     Creates a :dense NMatrix with a sequence of integers starting at
  #     zero until the matrix is filled. The parameters to the method
  #     are the dimensions of the matrix. Optionaly, one can specify a
  #     dtype as the last parameter (default is :float64).
  #     
  # Examples:
  # 
  #   seq(4) # =>   0.0   1.0   2.0   3.0
  #
  #   seq(3,3, :int32) # =>  0  1  2
  #                          3  4  5
  #                          6  7  8
  #

  def self.seq(*params)
  	dtype = params.last.is_a?(Symbol) ? params.pop   : :float64
	
    product = params.reduce(1) { |prod, n| prod *= n }
  	
  	NMatrix.new(params,  (0..prod-1).to_a, dtype )
  end

  # rand()
  #      
  #     Creates a :dense NMatrix with random numbers between 0 and 1 generated
  #     by Random::rand. The parameters to the method are the dimensions of the
  #     matrix.
  #     
  #     by Carlos Agarie <carlos@onox.com.br>
  #
  # Examples:
  # 
  #   rand(3,3) # => 0.4859439730644226   0.1783195585012436
  #                  0.23193766176700592  0.4503345191478729
  #

  def self.rand(*params)
    rng = Random.new
  
    product = params.reduce(1) { |prod, n| prod *= n }
  
    random_values = []
    product.times { |i| random_values << rng.rand }
  
    NMatrix.new params, random_values, :float32
  end

  ######################################
  #
  # FUNCIONS FOR MATLAB AND IDL REFUGEES
  #
  ######################################
  #
  # These are functions that replicate existing functionality, but
  # would probably be appreciated by MATLAB or IDL users.
  #

  # indgen() , findgen() , bindgen() , cindgen()
  #      
  #      These IDL functions are similar to seq() but less flexile.
  #      They produce on-dimensional matrices:
  #      
  #      indgen    -- Integer vector   --   seq(n, :int32)
  #      findgen   -- Float vector     --   seq(n, :float32)
  #      bindgen   -- Byte vector      --   seq(n, :byte)
  #      cindgen   -- Complex vector   --   seq(n, :complex64)
  #      

  def self.indgen(n)
  	seq(n, :int32)
  end

  def self.findgen(n)
  	seq(n, :float32)
  end

  def self.bindgen(n)
  	seq(n, :byte)
  end

  def self.cindgen(n)
  	seq(n, :complex64)
  end

  # linspace()
  #      
  #      This MATLAB function somewhat resembles seq(), but it differs
  #      enough that is likely to be legitimately useful to non-MATLAB
  #      refugees. This function takes three parameter, the last one
  #      being an integer.
  #      
  #      linspace( a, b, n )
  #      
  #      This returns a vector with n values from a to b.
  #      
  # Ex:  x = linspace(0, pi, 1000)
  #      y = sin(x)
  #

  def self.linspace(a, b, nsteps)
  	#
  	# Algorithm:  seq(n) * (b-a)/(n-1) + a
  	#
  	step = (b-a) * 1.0 / (nsteps - 1)
  	seq(nsteps) * NMatrix.new([nsteps], step) + NMatrix.new([nsteps], a)
  end
  
end