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
# == io/matlab/harwell_boeing.rb
#
# Harwell Boeing file reader (and eventually writer too).
# => Supports only assembled, non-symmetric, real matrices
# => Data types supported are exponential, floating point and integer
# => Returned NMatrix is of type :float64
#++

require_relative './fortran_format.rb'

class NMatrix
  module IO
    module HarwellBoeing

      class << self
        # Loads the contents of a valid Harwell Boeing format file and 
        # returns an NMatrix object with the values of the file and optionally
        # only the header info.
        # 
        # Supports only assembled, non-symmetric, real matrices. File name must
        # have matrix type as extension.
        # 
        # Example - test_file.rua
        # 
        # == Arguments
        # 
        # * +file_path+ - Path of the Harwell Boeing file  to load.
        # * +opts+      - Options for specifying whether you want
        #                 the values and  header or only the header.
        # 
        # == Options
        # 
        # * +:header+ - If specified as *true*, will return only the header of
        #               the HB file.Will return the NMatrix object and
        #               header as an array if left blank.
        # 
        # == Usage
        # 
        #   mat, head = NMatrix::IO::HarwellBoeing.load("test_file.rua")
        # 
        #   head = NMatrix::IO::HarwellBoeing.load("test_file.rua", {header: true})
        # 
        # == Alternate Usage
        # 
        # You can specify the file using NMatrix::IO::Reader.new("path/to/file")
        # and then call *header* or *values* on the resulting object.
        def load file_path, opts={}
          hb_obj = NMatrix::IO::HarwellBoeing::Reader.new(file_path)

          return hb_obj.header if opts[:header]

          [hb_obj.values, hb_obj.header]
        end
      end

      class Reader
        def initialize file_name
          raise(IOError, "Unsupported file format. Specify file as \
            file_name.rua.") if !file_name.match(/.*\.[rR][uU][aA]/)

          @file_name   = file_name
          @header      = {}
          @body        = nil
        end

        def header
          return @header if !@header.empty?
          @file = File.open @file_name, "r"

          line = @file.gets

          @header[:title] = line[0...72].strip
          @header[:key]   = line[72...80].strip

          line = @file.gets

          @header[:totcrd] = line[0...14] .strip.to_i
          @header[:ptrcrd] = line[14...28].strip.to_i
          @header[:indcrd] = line[28...42].strip.to_i
          @header[:valcrd] = line[42...56].strip.to_i
          @header[:rhscrd] = line[56...70].strip.to_i

          raise(IOError, "Right hand sides not supported.") \
           if @header[:rhscrd] > 0

          line = @file.gets

          @header[:mxtype] = line[0...3]

          raise(IOError, "Currently supports only real, assembled, unsymmetric \
            matrices.") if !@header[:mxtype].match(/RUA/)

          @header[:nrow]   = line[13...28].strip.to_i
          @header[:ncol]   = line[28...42].strip.to_i
          @header[:nnzero] = line[42...56].strip.to_i
          @header[:neltvl] = line[56...70].strip.to_i

          line = @file.gets

          fortran_reader = NMatrix::IO::FortranFormat::Reader

          @header[:ptrfmt] = fortran_reader.new(line[0...16].strip) .parse
          @header[:indfmt] = fortran_reader.new(line[16...32].strip).parse
          @header[:valfmt] = fortran_reader.new(line[32...52].strip).parse
          @header[:rhsfmt] = fortran_reader.new(line[52...72].strip).parse

          @header
        end

        def values
          @header      = header if @header.empty?
          @file.lineno = 5      if @file.lineno != 5
          @matrix      = NMatrix.new([ @header[:nrow], @header[:ncol] ], 
                                      0, dtype: :float64)

          read_column_pointers
          read_row_indices
          read_values

          @file.close
          
          assemble_matrix

          @matrix
        end

       private

        def read_column_pointers
          @col_ptrs  = []
          pointer_lines     = @header[:ptrcrd]
          pointers_per_line = @header[:ptrfmt][:repeat]
          pointer_width     = @header[:ptrfmt][:field_width]

          @col_ptrs = read_numbers :to_i, pointer_lines, pointers_per_line, 
                                             pointer_width

          @col_ptrs.map! {|c| c -= 1}
        end

        def read_row_indices
          @row_indices     = []
          row_lines        = @header[:indcrd]
          indices_per_line = @header[:indfmt][:repeat]
          row_width        = @header[:indfmt][:field_width]

          @row_indices = read_numbers :to_i, row_lines, indices_per_line, 
                                      row_width

          @row_indices.map! {|r| r -= 1}
        end

        def read_values
          @vals = []
          value_lines = @header[:valcrd]
          values_per_line = @header[:valfmt][:repeat]
          value_width    = @header[:valfmt][:field_width]

          @vals = read_numbers :to_f, value_lines, values_per_line, 
                                  value_width
        end

        def read_numbers to_dtype, num_of_lines, numbers_per_line, number_width
          data = []

          num_of_lines.times do 
            line  = @file.gets
            index = 0

            numbers_per_line.times do
              delimiter = index + number_width

              data << line[index...delimiter].strip.send(to_dtype)

              break if line.length <= delimiter
              index += number_width
            end
          end

          data
        end

        def assemble_matrix
          col = 0
          @col_ptrs[0..-2].each_index do |index|
            @col_ptrs[index].upto(@col_ptrs[index+1] - 1) do |row_ptr|
              row               = @row_indices[row_ptr]
              @matrix[row, col] = @vals[row_ptr]
            end

            col += 1
          end
        end
      end

    end
  end
end