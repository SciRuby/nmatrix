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
# == io/matlab/fortran_format.rb
#
# A parser for making sense of FORTRAN formats.
# => Only handles R (real), F (float) and E (exponential) format codes. 
#++

class NMatrix
  module IO
    module FortranFormat

      # Class for reading strings in FORTRAN format for specifying attributes
      # of numerical data in a file. Supports F (float), E (exponential) and 
      # R (real).
      # 
      # == Usage
      # 
      #   p = NMatrix::IO::FortranFormat::Reader.new("(16I5)")
      #   v = p.parse
      #   puts v #=> { :format_code => "INT_ID", 
      #          #=>   :repeat      =>       16,
      #          #=>   :field_width =>        5 }
      class Reader

        # Accepts a string in FORTRAN format and initializes the 
        # NMatrix::IO::FortranFormat::Reader object for further parsing of the 
        # data.
        # 
        # == Arguments
        # 
        # * +string+ - FORTRAN format string to be parsed.
        def initialize string
          @string = string
        end

        # Parses the FORTRAN format string passed in initialize and returns
        # a hash of the results.
        # 
        # == Result Hash Format
        # 
        # Take note that some of the below parameters may be absent in the hash
        # depending on the type of string being parsed.
        # 
        # * +:format_code+ - A string containing the format code of the read data. 
        #                    Can be "INT_ID", "FP_ID" or "EXP_ID" 
        # * +:repeat+      - Number of times this format will repeat in a line.
        # * +:field_width+ - Width of the numerical part of the number.
        # * +:post_decimal_width+ - Width of the numerals after the decimal point.
        # * +:exponent_width+ - Width of exponent part of the number.
        def parse
          raise(IOError, "Left or right parentheses missing") \
           if parentheses_missing? # change tests to handle 'raise' not return

          @result = {}
          @string = @string[1..-2]

          if valid_fortran_format?
            load_result
          else
            raise(IOError, "Invalid FORTRAN format specified. Only Integer, Float or Exponential acceptable.")
          end

          @result
        end

       private
        def parentheses_missing?
          true if @string[0] != '(' or @string[-1] != ')'
        end

        # Changing any of the following regular expressions can lead to disaster
        def valid_fortran_format?
          @mdata = @string.match(/\A(\d*)(I)(\d+)\z/) # check for integer format
          @mdata = @string.match(/\A(\d*)(F)(\d+)\.(\d+)\z/) \
           if @mdata.nil? # check for floating point if not integer
          @mdata =  @string.match(/\A(\d*)(E)(\d+)\.(\d+)(E)?(\d*)\z/) \
           if @mdata.nil? # check for exponential format if not floating point

          @mdata
        end

        def load_result
          if @mdata.to_a.include? "I"
            create_integer_hash
          elsif @mdata.to_a.include? "F"
            create_float_hash
          else
            create_exp_hash
          end
        end

        def create_integer_hash
          @result[:format_code] = "INT_ID"
          @result[:repeat]      = @mdata[1].to_i if !@mdata[1].empty?
          @result[:field_width] = @mdata[3].to_i
        end

        def create_float_hash
          @result[:format_code]        = "FP_ID"
          @result[:repeat]             = @mdata[1].to_i if !@mdata[1].empty?
          @result[:field_width]        = @mdata[3].to_i
          @result[:post_decimal_width] = @mdata[4].to_i
        end

        def create_exp_hash
          @result[:format_code]        = "EXP_ID"
          @result[:repeat]             = @mdata[1].to_i if !@mdata[1].empty?
          @result[:field_width]        = @mdata[3].to_i
          @result[:post_decimal_width] = @mdata[4].to_i
          @result[:exponent_width]     = @mdata[6].to_i if !@mdata[6].empty?
        end
      end
      
    end
  end
end