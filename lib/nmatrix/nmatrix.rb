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
# == nmatrix.rb
#
# This file adds a few additional pieces of functionality (e.g., inspect,
# pretty_print).

############
# Requires #
############

#######################
# Classes and Modules #
#######################

class NMatrix
	# Read and write extensions for NMatrix. These are only loaded when needed.
	module IO
		autoload :MatReader, 'nmatrix/io/mat_reader'
		autoload :Mat5Reader, 'nmatrix/io/mat5_reader'
	end

	# TODO: Make this actually pretty.
	def pretty_print(q = nil)
		raise NotImplementedError, 'Can only print rank 2 matrices.' unless rank == 2

    arr = []

		(0...shape[0]).each do |i|
			arr << (0...shape[1]).inject(Array.new) { |a, j| a << if o = self[i, j] then o else 'nil' end }
    end

    if q
      q.group(1, "", "\n") do
        q.seplist(arr, lambda { q.text "  " }, :each)  { |v| q.text v.to_s }
      end
    else
      puts arr.join("  ")
    end

	end
	alias :pp :pretty_print


	# Get the complex conjugate of this matrix. See also complex_conjugate! for
	# an in-place operation (provided the dtype is already :complex64 or
	# :complex128).
	#
	# Does not work on list matrices, but you can optionally pass in the type you
	# want to cast to if you're dealing with a list matrix.
	def complex_conjugate(new_stype = self.stype)
		self.cast(new_stype, NMatrix::upcast(dtype, :complex64)).complex_conjugate!
	end

	# Calculate the conjugate transpose of a matrix. If your dtype is already
	# complex, this should only require one copy (for the transpose).
	def conjugate_transpose
		self.transpose.complex_conjugate!
	end

	def hermitian?
		return false if self.rank != 2 or self.shape[0] != self.shape[1]
		
		if [:complex64, :complex128].include?(self.dtype)
			# TODO: Write much faster Hermitian test in C
			self.eql?(self.conjugate_transpose)
		else
			symmetric?
		end
	end

	def inspect
		original_inspect = super()
		original_inspect = original_inspect[0...original_inspect.size-1]
		original_inspect + inspect_helper.join(" ") + ">"
	end

	def __yale_ary__to_s(sym)
		ary = self.send("__yale_#{sym.to_s}__".to_sym)
		
		'[' + ary.collect { |a| a ? a : 'nil'}.join(',') + ']'
	end

	class << self
		def cblas_gemm(a, b, c = nil, alpha = 1.0, beta = 0.0, transpose_a = false, transpose_b = false, m = nil, n = nil, k = nil, lda = nil, ldb = nil, ldc = nil)
			raise ArgumentError, 'Expected dense NMatrices as first two arguments.' unless a.is_a?(NMatrix) and b.is_a?(NMatrix) and a.stype == :dense and b.stype == :dense
			raise ArgumentError, 'Expected nil or dense NMatrix as third argument.' unless c.nil? or (c.is_a?(NMatrix) and c.stype == :dense)
			raise ArgumentError, 'NMatrix dtype mismatch.'													unless a.dtype == b.dtype and (c ? a.dtype == c.dtype : true)

			# First, set m, n, and k, which depend on whether we're taking the
			# transpose of a and b.
			if c
				m ||= c.shape[0]
				n ||= c.shape[1]
				k ||= transpose_a ? a.shape[0] : a.shape[1]
				
			else
				if transpose_a
					# Either :transpose or :complex_conjugate.
					m ||= a.shape[1]
					k ||= a.shape[0]
					
				else
					# No transpose.
					m ||= a.shape[0]
					k ||= a.shape[1]
				end
				
				n ||= transpose_b ? b.shape[0] : b.shape[1]
				c		= NMatrix.new([m, n], a.dtype)
			end

			# I think these are independent of whether or not a transpose occurs.
			lda ||= a.shape[1]
			ldb ||= b.shape[1]
			ldc ||= c.shape[1]

			# NM_COMPLEX64 and NM_COMPLEX128 both require complex alpha and beta.
			if a.dtype == :complex64 or a.dtype == :complex128
				alpha = Complex.new(1.0, 0.0) if alpha == 1.0
				beta  = Complex.new(0.0, 0.0) if beta  == 0.0
			end

			# For argument descriptions, see: http://www.netlib.org/blas/dgemm.f
			NMatrix.__cblas_gemm__(transpose_a, transpose_b, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

			return c
		end

		def cblas_gemv(a, x, y = nil, alpha = 1.0, beta = 0.0, transpose_a = false, m = nil, n = nil, lda = nil, incx = nil, incy = nil)
			m ||= transpose_a ? a.shape[1] : a.shape[0]
			n ||= transpose_a ? a.shape[0] : a.shape[1]
			
			lda		||= a.shape[1]
			incx	||= 1
			incy	||= 1
			
			# NM_COMPLEX64 and NM_COMPLEX128 both require complex alpha and beta.
			if a.dtype == :complex64 or a.dtype == :complex128
				alpha = Complex.new(1.0, 0.0) if alpha == 1.0
				beta  = Complex.new(0.0, 0.0) if beta  == 0.0
			end
			
			NMatrix.__cblas_gemv__(transpose_a, m, n, alpha, a, lda, x, incx, beta, y, incy)
			
			return y
		end
		
		def load_file(file_path)
			NMatrix::IO::Mat5Reader.new(File.open(file_path, 'rb')).to_ruby
		end
		
		# Helper function for loading a file in the first sparse format given here:
		#   http://math.nist.gov/MatrixMarket/formats.html
		#
		# Override type specifier (e.g., 'real') using :read_with => :to_f (or any other string-to-numeric conversion
		# function), and with :dtype => :float32 or :dtype => :int8 to force storage in a lesser type.
		def load_matrix_matrix_coordinate_file(filename, options = {})
			f = File.new(filename, "r")

			func	= options[:read_with]
			dtype = options[:dtype]
			
			line = f.gets
			raise IOError, 'Incorrect file type specifier.' unless line =~ /^%%MatrixMarket\ matrix\ coordinate/
			spec = line.split
			
			case spec[3]
			when 'real'
				func	||= :to_f
				dtype ||= :float64
			when 'integer'
				func	||= :to_i
				dtype ||= :int64
			when 'complex'
				func	||= :to_complex
				dtype ||= :complex128
			when 'rational'
				func	||= :to_rational
				dtype ||= :rational128
			else
				raise ArgumentError, 'Unrecognized dtype.'
			end unless func and dtype
			
			begin
				line = f.gets
			end while line =~ /^%/
			
			# Close the file.
			f.close
			
			rows, cols, entries = line.split.collect { |x| x.to_i }
			
			matrix = NMatrix.new(:yale, [rows, cols], entries, dtype)
			
			entries.times do
				i, j, v = line.split
				matrix[i.to_i - 1, j.to_i - 1] = v.send(func)
			end
			
			matrix
		end
	end

	protected
	def inspect_helper
		ary = []
		ary << "shape:[#{shape.join(',')}]" << "dtype:#{dtype}" << "stype:#{stype}"

		if stype == :yale
			ary <<	"capacity:#{capacity}" << "ija:#{__yale_ary__to_s(:ija)}" << "ia:#{__yale_ary__to_s(:ia)}" <<
							"ja:#{__yale_ary__to_s(:ja)}" << "a:#{__yale_ary__to_s(:a)}" << "d:#{__yale_ary__to_s(:d)}" <<
							"lu:#{__yale_ary__to_s(:lu)}" << "yale_size:#{__yale_size__}"
		end

		ary
	end
end

#######################################
#
# SIMPLE IN-LINE CONSTRUCTOR
#
# by Daniel Carrera <dcarrera@hush.com>
#
#######################################
# 
# NMatrix needs to have a succinct way to create a matrix by specifying
# the components directly. This is very usefeul for using NMatrix as an
# advanced calculator, it is useful for learning NMatrix and it is also
# useful for testing language features or developing algorithms.
# 
# The N[] function provides a way to create a matrix in a way that is
# very short and very natural, simply by specifying the components in
# the traditional Ruby array syntax.  Optionaly, one can specify a
# dtype as the last parameter (default is :float64).
# 
# a = N[ 1,2,3,4 ]          =>  1.0  2.0  3.0  4.0
# 
# a = N[ 1,2,3,4, :int32 ]  =>  1  2  3  4
# 
# a = N[ [1,2,3], [3,4,5] ] =>  1.0  2.0  3.0
#                               3.0  4.0  5.0
# 
# 
# SYNTAX COMPARISON:
# 
#     MATLAB:		a = [ [1 2 3] ; [4 5 6] ]   or  [ 1 2 3 ; 4 5 6 ]
#     IDL:			a = [ [1,2,3] , [4,5,6] ]
#     NumPy:		a = array( [1,2,3], [4,5,6] )
#     
#     SciRuby:      a = N[ [1,2,3], [4,5,6] ]
#     Ruby array:   a =  [ [1,2,3], [4,5,6] ]
#

class N

	def N.[](*params)
		dtype = params.last.is_a?(Symbol) ? params.pop : :float64
		
		#
		# First find the dimensions of the array.
		#
		i = 0
		dim = []
		foo = params
		while foo.is_a?(Array)
			dim[i] = foo.length
			foo    = foo[0]
			i += 1
		end
		
		#
		# Then flatten the array.
		#
		NMatrix.new( dim, params.flatten, dtype )
	end
end