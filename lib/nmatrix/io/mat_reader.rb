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
# == io/mat_reader.rb
#
# Base class for .mat file reading (Matlab files).
#
#++

require 'packable'

module NMatrix::IO::Matlab

  # Class for parsing a .mat file stream.
  #
  # The full format of .mat files is available here:
  # * http://www.mathworks.com/help/pdf_doc/matlab/matfile_format.pdf
  class MatReader #:nodoc:
    MDTYPE_UNPACK_ARGS = {
      :miINT8   => [Integer, {:signed    => true,    :bytes => 1}],
      :miUINT8  => [Integer, {:signed    => false,   :bytes => 1}],
      :miINT16  => [Integer, {:signed    => true,    :bytes => 2}],
      :miUINT16 => [Integer, {:signed    => false,   :bytes => 2}],
      :miINT32  => [Integer, {:signed    => true,    :bytes => 4}],
      :miUINT32 => [Integer, {:signed    => false,   :bytes => 4}],
      :miSINGLE => [Float,   {:precision => :single,
                              :bytes => 4, :endian => :native}],
      :miDOUBLE => [Float,   {:precision => :double,
                              :bytes => 4, :endian => :native}],
      :miINT64  => [Integer, {:signed    => true,    :bytes => 8}],
      :miUINT64 => [Integer, {:signed    => false,   :bytes => 8}]
    }

    DTYPE_PACK_ARGS = {
      :byte       => [Integer, {:signed => false,
                                :bytes => 1}],
      :int8       => [Integer, {:signed => true,
                                :bytes => 1}],
      :int16      => [Integer, {:signed => true,
                                :bytes => 2}],
      :int32      => [Integer, {:signed => true,
                                :bytes => 4}],
      :int64      => [Integer, {:signed => true,
                                :bytes => 8}],
      :float32    => [Float,   {:precision => :single,
                                :bytes => 4, :endian => :native}],
      :float64    => [Float,   {:precision => :double,
                                :bytes => 8, :endian => :native}],
      :complex64  => [Float,   {:precision => :single,
                                :bytes => 4, :endian => :native}], #2x
      :complex128 => [Float,   {:precision => :double,
                                :bytes => 8, :endian => :native}]
    }

    ITYPE_PACK_ARGS = {
      :uint8  => [Integer, {:signed => false, :bytes => 1}],
      :uint16 => [Integer, {:signed => false, :bytes => 2}],
      :uint32 => [Integer, {:signed => false, :bytes => 4}],
      :uint64 => [Integer, {:signed => false, :bytes => 8}],
    }

    NO_REPACK = [:miINT8, :miUINT8, :miINT16,
                 :miINT32, :miSINGLE, :miDOUBLE, :miINT64]

    # Convert from MATLAB dtype to NMatrix dtype.
    MDTYPE_TO_DTYPE = {
      :miUINT8  => :byte,
      :miINT8   => :int8,
      :miINT16  => :int16,
      :miUINT16 => :int16,
      :miINT32  => :int32,
      :miUINT32 => :int32,
      :miINT64  => :int64,
      :miUINT64 => :int64,
      :miSINGLE => :float32,
      :miDOUBLE => :float64
    }

    MDTYPE_TO_ITYPE = {
      :miUINT8  => :uint8,
      :miINT8   => :uint8,
      :miINT16  => :uint16,
      :miUINT16 => :uint16,
      :miINT32  => :uint32,
      :miUINT32 => :uint32,
      :miINT64  => :uint64,
      :miUINT64 => :uint64
    }

    # Before release v7.1 (release 14) matlab (TM) used the system
    # default character encoding scheme padded out to 16-bits. Release 14
    # and later use Unicode. When saving character data, R14 checks if it
    # can be encoded in 7-bit ascii, and saves in that format if so.
    MDTYPES = [
               nil,
               :miINT8,
               :miUINT8,
               :miINT16,
               :miUINT16,
               :miINT32,
               :miUINT32,
               :miSINGLE,
               nil,
               :miDOUBLE,
               nil,
               nil,
               :miINT64,
               :miUINT64,
               :miMATRIX,
               :miCOMPRESSED,
               :miUTF8,
               :miUTF16,
               :miUTF32
              ]

    MCLASSES = [
                nil,
                :mxCELL,
                :mxSTRUCT,
                :mxOBJECT,
                :mxCHAR,
                :mxSPARSE,
                :mxDOUBLE,
                :mxSINGLE,
                :mxINT8,
                :mxUINT8,
                :mxINT16,
                :mxUINT16,
                :mxINT32,
                :mxUINT32,
                :mxINT64,
                :mxUINT64,
                :mxFUNCTION,
                :mxOPAQUE,
                :mxOBJECT_CLASS_FROM_MATRIX_H
               ]

    attr_reader :byte_order


    # call-seq:
    #     new(stream, options = {}) -> MatReader
    #
    # * *Raises* :
    #   - +ArgumentError+ -> First argument must be IO.
    #
    def initialize(stream, options = {})
      raise ArgumentError, 'First arg must be IO.' unless stream.is_a?(::IO)

      @stream     = stream
      @byte_order = options[:byte_order] || guess_byte_order
    end

    # call-seq:
    #     guess_byte_order -> Symbol
    #
    def guess_byte_order
      # Assume native, since we don't know what type of file we have.
      :native
    end

    protected

    attr_reader :stream
  end
end
