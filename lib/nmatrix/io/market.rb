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
# SciRuby is Copyright (c) 2010 - 2016, Ruby Science Foundation
# NMatrix is Copyright (c) 2012 - 2016, John Woods and the Ruby Science Foundation
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
# == io/market.rb
#
# MatrixMarket reader and writer.
#
#++

# Matrix Market is a repository of test data for use in studies of algorithms
# for numerical linear algebra. There are 3 file formats used:
#
# - Matrix Market Exchange Format.
# - Harwell-Boeing Exchange Format.
# - Coordinate Text File Format. (to be phased out)
#
# This module can load and save the first format. We might support
# Harwell-Boeing in the future.
#
# The MatrixMarket format is documented in:
# * http://math.nist.gov/MatrixMarket/formats.html
module NMatrix::IO::Market
  CONVERTER_AND_DTYPE = {
    :real => [:to_f, :float64],
    :complex => [:to_c, :complex128],
    :integer => [:to_i, :int64],
    :pattern => [:to_i, :byte]
  } #:nodoc:

  ENTRY_TYPE = {
    :byte => :integer, :int8 => :integer, :int16 => :integer,
    :int32 => :integer, :int64 => :integer,:float32 => :real,
    :float64 => :real, :complex64 => :complex, :complex128 => :complex
  } #:nodoc:

  class << self

    # call-seq:
    #     load(filename) -> NMatrix
    #
    # Load a MatrixMarket file. Requires a +filename+ as an argument.
    #
    # * *Arguments* :
    #   - +filename+ -> String with the filename to be saved.
    # * *Raises* :
    #   - +IOError+ -> expected type code line beginning with '%%MatrixMarket matrix'
    def load(filename)

      f = File.new(filename, "r")

      header = f.gets
      header.chomp!
      raise(IOError, "expected type code line beginning with '%%MatrixMarket matrix'") \
       if header !~ /^\%\%MatrixMarket\ matrix/

      header = header.split

      entry_type = header[3].downcase.to_sym
      symmetry   = header[4].downcase.to_sym
      converter, default_dtype = CONVERTER_AND_DTYPE[entry_type]

      if header[2] == 'coordinate'
        load_coordinate f, converter, default_dtype, entry_type, symmetry
      else
        load_array f, converter, default_dtype, entry_type, symmetry
      end
    end

    # call-seq:
    #     save(matrix, filename, options = {}) -> true
    #
    # Can optionally set :symmetry to :general, :symmetric, :hermitian; and can
    # set :pattern => true if you're writing a sparse matrix and don't want
    # values stored.
    #
    # * *Arguments* :
    #   - +matrix+ -> NMatrix with the data to be saved.
    #   - +filename+ -> String with the filename to be saved.
    # * *Raises* :
    #   - +DataTypeError+ -> MatrixMarket does not support Ruby objects.
    #   - +ArgumentError+ -> Expected two-dimensional NMatrix.
    def save(matrix, filename, options = {})
      options = {:pattern => false,
        :symmetry => :general}.merge(options)

      mode = matrix.stype == :dense ? :array : :coordinate
      if [:object].include?(matrix.dtype)
        raise(DataTypeError, "MatrixMarket does not support Ruby objects")
      end
      entry_type = options[:pattern] ? :pattern : ENTRY_TYPE[matrix.dtype]

      raise(ArgumentError, "expected two-dimensional NMatrix") \
       if matrix.dim != 2

      f = File.new(filename, 'w')

      f.puts "%%MatrixMarket matrix #{mode} #{entry_type} #{options[:symmetry]}"

      if matrix.stype == :dense
        save_array matrix, f, options[:symmetry]
      elsif [:list,:yale].include?(matrix.stype)
        save_coordinate matrix, f, options[:symmetry], options[:pattern]
      end

      f.close

      true
    end


    protected

    def save_coordinate matrix, file, symmetry, pattern
      # Convert to a hash in order to store
      rows = matrix.to_h

      # Count non-zeros
      count = 0
      rows.each_pair do |i, columns|
        columns.each_pair do |j, val|
          next if symmetry != :general && j > i
          count += 1
        end
      end

      # Print dimensions and non-zeros
      file.puts "#{matrix.shape[0]}\t#{matrix.shape[1]}\t#{count}"

      # Print coordinates
      rows.each_pair do |i, columns|
        columns.each_pair do |j, val|
          next if symmetry != :general && j > i
          file.puts(pattern ? "\t#{i+1}\t#{j+1}" : "\t#{i+1}\t#{j+1}\t#{val}")
        end
      end

      file
    end


    def save_array matrix, file, symmetry
      file.puts [matrix.shape[0], matrix.shape[1]].join("\t")

      if symmetry == :general
        (0...matrix.shape[1]).each do |j|
          (0...matrix.shape[0]).each do |i|
            file.puts matrix[i,j]
          end
        end
      else # :symmetric, :'skew-symmetric', :hermitian
        (0...matrix.shape[1]).each do |j|
          (j...matrix.shape[0]).each do |i|
            file.puts matrix[i,j]
          end
        end
      end

      file
    end


    def load_array file, converter, dtype, entry_type, symmetry
      mat = nil

      line = file.gets
      line.chomp!
      line.lstrip!

      fields = line.split

      mat = NMatrix.new :dense, [fields[0].to_i, fields[1].to_i], dtype

      (0...mat.shape[1]).each do |j|
        (0...mat.shape[0]).each do |i|
          datum = file.gets.chomp.send(converter)
          mat[i,j] = datum

          unless i == j || symmetry == :general
            if symmetry == :symmetric
              mat[j,i] = datum
            elsif symmetry == :hermitian
              mat[j,i] = Complex.new(datum.real, -datum.imag)
            elsif symmetry == :'skew-symmetric'
              mat[j,i] = -datum
            end
          end
        end
      end

      file.close

      mat
    end


    # Creates a :list NMatrix from a coordinate-list MatrixMarket file.
    def load_coordinate file, converter, dtype, entry_type, symmetry

      mat = nil

      # Read until we get the dimensions and nonzeros
      while line = file.gets
        line.chomp!
        line.lstrip!
        line, comment = line.split('%', 2) # ignore comments
        if line.size > 4
          shape0, shape1 = line.split
          mat = NMatrix.new(:list, [shape0.to_i, shape1.to_i], 0, dtype)
          break
        end
      end

      # Now read the coordinates
      while line = file.gets
        line.chomp!
        line.lstrip!
        line, comment = line.split('%', 2) # ignore comments

        next unless line.size >= 5 # ignore empty lines

        fields = line.split

        i = fields[0].to_i - 1
        j = fields[1].to_i - 1
        datum = entry_type == :pattern ? 1 : fields[2].send(converter)

        mat[i, j] = datum # add to the matrix
        unless i == j || symmetry == :general
          if symmetry == :symmetric
            mat[j, i] = datum
          elsif symmetry == :'skew-symmetric'
            mat[j, i] = -datum
          elsif symmetry == :hermitian
            mat[j, i] = Complex.new(datum.real, -datum.imag)
          end
        end
      end

      file.close

      mat
    end
  end
end
