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
module NMatrix::IO::Csv
  # Reader for a version .csv file.

  class << self
    # call-seq:
    #     load(filename) -> NMatrix
    #     load(filename, dtype) -> NMatrix
    #     load(filename, delim) -> NMatrix
    #     load(filename, dtype, delim) -> NMatrix
    #
    # * *Arguments* :
    #   - +filename+ -> String giving the name of the file to be loaded.
    #   - +dtype+ -> Data type expected for the data of the file, 
		#   :float64 as default
    #   - +delim+ -> Delimiter of the csv file, comma as default
    #
    # Load a Comma Seperated Value CSV file as a matrix.
    def load(filename, *opts)
      CsvReader.new(filename, *opts).matrix
    end
  end

  class CsvReader #:nodoc:
    @@csv_plugged = true
    begin
      require_relative '../csv.rb'
    rescue
      @@csv_plugged = false
    end

    # call-seq:
    #     Csv::CsvReader.new(filename) -> CsvReader
    #
    # * *Arguments* :
    #   - +filename+ -> String giving the name of the file to be loaded.
    # * *Raises* :
    #   - +IOError+ -> CSV file not Exist
    #
    def initialize filename, *opts

      raise(IOError, "File does not exist.") unless File.exist?(filename)

      @matrix = read_matrix filename, *opts # Ensure it is a hash
    end

    attr_reader :matrix
    protected
      # Dispatch the read action
      # to actual read functions 
      def read_matrix filename, *opts
        params = opts[-1]
        params = {} if !params
        delim = params[:delim] || ','
        target_type = params[:dtype] || :object
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
			# Maybe something wrong will happen if the fields are quoted
      def ruby_read filename, delim, target_type
        rows = []
        result_arr = []
        dim2 = 0;
        File.new(filename).readlines.each do |line|
          rows << line.split(delim).map! do |str| str.strip! end
          dim2 = [dim2, rows[-1].size].max
        end
        rows.map! do |row|
          row.concat(Array.new(dim2 - row.size, nil))
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
