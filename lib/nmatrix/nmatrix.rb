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

require_relative './shortcuts.rb'

#######################
# Classes and Modules #
#######################

class NMatrix
	# Read and write extensions for NMatrix. These are only loaded when needed.
	module IO
    module Matlab
      class << self
        def load_mat file_path
          NMatrix::IO::Matlab::Mat5Reader.new(File.open(file_path, "rb+")).to_ruby
        end
      end

      # FIXME: Remove autoloads
      autoload :MatReader, 'nmatrix/io/mat_reader'
      autoload :Mat5Reader, 'nmatrix/io/mat5_reader'
    end
	end

	# TODO: Make this actually pretty.
	def pretty_print(q = nil)
		if dim != 2 || (dim == 2 && shape[1] > 10) # FIXME: Come up with a better way of restricting the display
      inspect
    else

      arr = (0...shape[0]).map do |i|
        ary = []
        (0...shape[1]).each do |j|
          o = begin
            self[i, j]
          rescue ArgumentError
            nil
          end
          ary << (o.nil? ? 'nil' : o)
        end
        ary.inspect
      end

      if q.nil?
        puts arr.join("\n")
      else
        q.group(1, "", "\n") do
          q.seplist(arr, lambda { q.text "  " }, :each)  { |v| q.text v.to_s }
        end
      end

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
		return false if self.dim != 2 or self.shape[0] != self.shape[1]
		
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
			ary <<	"capacity:#{capacity}"

      # These are enabled by the DEBUG_YALE compiler flag in extconf.rb.
      if respond_to?(:__yale_a__)
        ary << "ija:#{__yale_ary__to_s(:ija)}" << "ia:#{__yale_ary__to_s(:ia)}" <<
				  			"ja:#{__yale_ary__to_s(:ja)}" << "a:#{__yale_ary__to_s(:a)}" << "d:#{__yale_ary__to_s(:d)}" <<
					  		"lu:#{__yale_ary__to_s(:lu)}" << "yale_size:#{__yale_size__}"
      end

		end

		ary
	end
end