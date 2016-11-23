#--
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
# SciRuby is Copyright (c) 2010 - 2014, Ruby Science Foundation
# NMatrix is Copyright (c) 2012 - 2014, John Woods and the Ruby Science Foundation
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
# == io/csv_reader.rb
#
# a csv file reader.
#
#++

# TODO: ADD quote sepration
module NMatrix::IO::Csv
  # Reader for a version .csv file.

  class << self
    # call-seq:
    #     load(filename) -> NMatrix
    #     load(filename, options) -> NMatrix
    #
    # Load a Comma Seperated Value CSV file as a matrix.
    #
    # * *Arguments* :
    #   - +filename+ -> String giving the name of the file to be loaded.
    #   - +options+ -> Hash with options (:delim String, ',' by default; :dtype Symbol, 
    #   may be the type supported by SciRuby, :float64 by default)
    #
    # * *Returns* :
    #   - A NMatrix::IO::Csv::CsvReader object 
    #   with the matrix stored in its +matrix+ property
    #
    # * *Using* *options* :
    #   load('mat.csv', delim: ',')
    #   load('mat.csv', delim: ',' , dtype: :int32)
    #
    #   PLEASE NOTE THAT: 
    #   specifying :dtype with :object will result in a :float64 cast
    #
    def load(filename, options = {})
      CsvReader.new(filename, options).matrix
    end
  end

  class CsvReader #:nodoc:
    @@csv_plugged = true
    if File.exist?('lib/nmatrix/csv.rb')
      require_relative '../csv.rb'
    else
      @@csv_plugged = false
      require 'csv.rb' # csv parser from stdlib
    end

    # call-seq:
    #     Csv::CsvReader.new(filename) -> CsvReader
    #
    # * *Arguments* :
    #   - +filename+ -> String giving the name of the file to be loaded.
    #   - +options+ -> Hash with options (:delim String, ',' by default; :dtype Symbol, 
    #   may be the type supported by SciRuby, :float64 by default)
    #
    # * *Returns* :
    #   - A NMatrix::IO::Csv::CsvReader object 
    #   with the matrix stored in its +matrix+ property
    #
    # * *Raises* :
    #   - +IOError+ -> CSV file not Exist
    #
    #   PLEASE NOTE THAT: 
    #   specifying :dtype with :object will result in a :float64 cast
    #
    def initialize filename, options = {}

      raise(IOError, "File does not exist.") unless File.exist?(filename)

      @matrix = read_matrix filename, options
    end

    attr_reader :matrix
    protected
      # Dispatch the read action
      # to actual read functions 
      def read_matrix filename, options = {}
        delim = options[:delim] || ','
        target_type = options[:dtype] || :object
        if @@csv_plugged   # libcsv used
          libcsv_read filename, delim, target_type
        else # pure ruby reader
          ruby_read filename, delim, target_type
        end
      end

      # Use libcsv to read the matrix
      def libcsv_read filename, delim, target_type
        result = Csv.parse_file filename, delim
        if !result
          @matrix = nil
        else
          cast_raw_data! result[1], target_type
          @matrix = NMatrix.new result[0], result[1], dtype: target_type
        end
      end

      # Use pure ruby to read the matrix
      # Using CSV from ruby stdlib
      def ruby_read filename, delim, target_type
        result_arr = []
        rows = CSV.read(filename, col_sep: delim)
        dim2 = rows.map do |row| row.size end.max
        rows.map! do |row|
          row.concat(Array.new(dim2 - row.size, nil)) # fill unaligned elements with nil
          result_arr.concat(row)
        end
        cast_raw_data! result_arr, target_type
        @matrix = NMatrix.new [rows.size, dim2], result_arr, dtype: target_type
      end

      # Cast elements of data into given target_type
      # if target_type is :object, then elements are converted into float
      def cast_raw_data! data, target_type
        case target_type
          when :int32, :int64
            data.map! do |v| v = v.to_i end
          when :float32, :float64, :object
            data.map! do |v| v = v.to_f end
          when :complex64, :complex128
            data.map! do |v| v = v.to_c end
          else # if nothing specified or object
            data.map! do |v|
              v = v.to_f
            end
            target_type = :float64
        end
      end
  end
end
