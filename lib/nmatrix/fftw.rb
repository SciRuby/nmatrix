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
# == fftw.rb
#
# ruby file for the nmatrix-fftw gem. Loads the C extension and defines
# nice ruby interfaces for FFTW functions.
#++

require 'nmatrix/nmatrix.rb'
require "nmatrix_fftw.so"

class NMatrix

  # Compute 1D FFT of the matrix using FFTW default parameters.
  # @return [NMatrix] NMatrix of dtype :complex128 containing computed values.
  # @example Compute 1D FFT of an NMatrix.
  #   nm = NMatrix.new([10],
  #     [
  #       Complex(9.32,0), Complex(44,0), Complex(125,0), Complex(34,0),
  #       Complex(31,0),   Complex(44,0), Complex(12,0),  Complex(1,0),
  #       Complex(53.23,0),Complex(-23.23,0)
  #     ], dtype: :complex128)
  #   nm.fft
  def fft
    input = self.dtype == :complex128 ? self : self.cast(dtype: :complex128)
    plan  = NMatrix::FFTW::Plan.new([self.size])
    plan.set_input input
    plan.execute
    plan.output
  end

  # Compute 2D FFT of a 2D matrix using FFTW default parameters.
  # @return [NMatrix] NMatrix of dtype :complex128 containing computed values.
  def fft2
    raise ShapeError, "Shape must be 2 (is #{self.shape})" if self.shape.size != 2
    input = self.dtype == :complex128 ? self : self.cast(dtype: :complex128)
    plan  = NMatrix::FFTW::Plan.new(self.shape, dim: 2)
    plan.set_input input
    plan.execute
    plan.output
  end

  module FFTW
    class Plan
      # Hash which holds the numerical values of constants that determine
      # the kind of transform that will be computed for a real input/real
      # output instance. These are one-one mappings to the respective constants
      # specified in FFTW. For example, for specifying the FFTW_R2HC constant
      # as the 'kind', pass the symbol :r2hc.
      #
      # @see http://www.fftw.org/fftw3_doc/Real_002dto_002dReal-Transform-Kinds.html#Real_002dto_002dReal-Transform-Kinds
      REAL_REAL_FFT_KINDS_HASH = {
        r2hc:    0,
        hc2r:    1,
        dht:     2,
        redft00: 3,
        redft01: 4,
        redft10: 5,
        redft11: 6,
        rodft00: 7,
        rodft01: 9,
        rodft10: 8,
        rodft11: 10
      }

      # Hash holding the numerical values of the flags that are passed in the 
      # `flags` argument of a FFTW planner routine. Multiple flags can be passed
      # to one instance of the planner. Their values are OR'd ('|') and then passed.
      # For example, for passing the FFTW_ESTIMATE constant, use :estimate.
      #
      # nmatrix-fftw supports the following flags into the planning routine:
      # * :estimate - Equivalent to FFTW_ESTIMATE. Specifies that, instead of 
      #   actual measurements of different algorithms, a simple heuristic is 
      #   used to pick a (probably sub-optimal) plan quickly. With this flag, 
      #   the input/output arrays are not overwritten during planning.
      # * :measure - Equivalent to FFTW_MEASURE. Tells FFTW to find an optimized
      #   plan by actually computing several FFTs and measuring their execution
      #   time. Depending on your machine, this can take some time (often a few 
      #   seconds).
      # * :patient - Equivalent to FFTW_PATIENT. Like FFTW_MEASURE, but considers
      #   a wider range of algorithms and often produces a “more optimal” plan 
      #   (especially for large transforms), but at the expense of several times
      #   longer planning time (especially for large transforms).
      # * :exhaustive - Equivalent to FFTW_EXHAUSTIVE. Like FFTW_PATIENT, but 
      #   considers an even wider range of algorithms, including many that we 
      #   think are unlikely to be fast, to produce the most optimal plan but 
      #   with a substantially increased planning time.
      #
      # @see http://www.fftw.org/fftw3_doc/Planner-Flags.html#Planner-Flags
      FLAG_VALUE_HASH = {
        estimate: 64,
        measure: 0,
        exhaustive: 8,
        patient: 32
      }

      # Hash holding numerical values of the direction in which a :complex_complex
      # type FFT should be performed.
      #
      # @see http://www.fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html#Complex-One_002dDimensional-DFTs
      # (The fourth argument, sign, can be either FFTW_FORWARD (-1) or 
      # FFTW_BACKWARD (+1), and indicates the direction of the transform you are
      # interested in; technically, it is the sign of the exponent in the transform)
      FFT_DIRECTION_HASH = {
        forward: -1,
        backward: 1
      }

      # Hash holding numerical equivalents of the DFT type. Used for determining
      # DFT type in C level.
      DATA_TYPE_HASH = {
        complex_complex: 0,
        real_complex:    1,
        complex_real:    2,
        real_real:       3
      }

      # Array holding valid options that can be passed into NMatrix::FFTW::Plan
      # so that invalid options aren't passed.
      VALID_OPTS = [:dim, :type, :direction, :flags, :real_real_kind]

      # @!attribute [r] shape
      #   @return [Array] Shape of the plan. Sequence of Fixnums.
      attr_reader :shape

      # @!attribute [r] size
      #   @return [Numeric] Size of the plan.
      attr_reader :size

      # @!attribute [r] type
      #   @return [Symbol] Type of the plan. Can be :complex_complex, 
      #   :complex_real, :real_complex or :real_real
      attr_reader :type

      # @!attribute [r] direction
      #   @return [Symbol] Can be :forward of :backward. Indicates the direction
      #   of the transform you are interested in; technically, it is the sign of
      #   the exponent in the transform. Valid only for :complex_complex type.
      attr_reader :direction

      # @!attribute [r] flags
      #   @return [Array<Symbol>] Can contain one or more symbols from
      #   FLAG_VALUE_HASH. Determines how the planner is prepared.
      #   @see FLAG_VALUE_HASH
      attr_reader :flags

      # @!attribute [r] dim
      #   @return [Fixnum] Dimension of the FFT. Should be 1 for 1-D FFT, 2 for
      #   2-D FFT and so on.
      attr_reader :dim

      # @!attribute [r] input
      #   @return [NMatrix] Input NMatrix. Will be valid once the 
      #   NMatrix::FFTW::Plan#set_input method has been called.
      attr_reader :input

      # @!attribute [r] output
      #   @return [NMatrix] Output NMatrix. Will be valid once the 
      #   NMatrix::FFTW::Plan#execute method has been called.
      attr_reader :output

      # @!attribute [r] real_real_kind
      #   @return [Symbol] Specifies the kind of real to real FFT being performed.
      #   This is a symbol from REAL_REAL_FFT_KINDS_HASH. Only valid when type
      #   of transform is of type :real_real.
      #   @see REAL_REAL_FFT_KINDS_HASH
      #   @see http://www.fftw.org/fftw3_doc/Real_002dto_002dReal-Transform-Kinds.html#Real_002dto_002dReal-Transform-Kinds
      attr_reader :real_real_kind

      # Create a plan for a DFT. The FFTW library requires that you first create
      # a plan for performing a DFT, so that FFTW can optimize its algorithms
      # according to your computer's hardware and various user supplied options.
      # 
      # @see http://www.fftw.org/doc/Using-Plans.html 
      #   For a comprehensive explanation of the FFTW planner.
      # @param shape [Array, Fixnum] Specify the shape of the plan. For 1D
      #   fourier transforms this can be a single number specifying the length of 
      #   the input. For multi-dimensional transforms, specify an Array containing
      #   the length of each dimension.
      # @param [Hash] opts the options to create a message with.
      # @option opts [Fixnum] :dim (1) The number of dimensions of the Fourier
      #   transform. If 'shape' has more numbers than :dim, the number of dimensions
      #   specified by :dim will be considered when making the plan.
      # @option opts [Symbol] :type (:complex_complex) The type of transform to
      #   perform based on the input and output data desired. The default value
      #   indicates that a transform is being planned that uses complex numbers
      #   as input and generates complex numbers as output. Similarly you can
      #   use :complex_real, :real_complex or :real_real to specify the kind
      #   of input and output that you will be supplying to the plan.
      #   @see DATA_TYPE_HASH
      # @option opts [Symbol, Array] :flags (:estimate) Specify one or more flags
      #   which denote the methodology that is used for deciding the algorithm used
      #   when planning the fourier transform. Use one or more of :estimate, :measure,
      #   :exhaustive and :patient. These flags map to the planner flags specified
      #   at http://www.fftw.org/fftw3_doc/Planner-Flags.html#Planner-Flags.
      #   @see REAL_REAL_FFT_KINDS_HASH
      # @option opts [Symbol] :direction (:forward) The direction of a DFT of
      #   type :complex_complex. Technically, it is the sign of the exponent in 
      #   the transform. :forward corresponds to -1 and :backward to +1.
      #   @see FFT_DIRECTION_HASH
      # @option opts [Array] :real_real_kind When the type of transform is :real_real,
      #   specify the kind of transform that should be performed FOR EACH AXIS
      #   of input. The position of the symbol in the Array corresponds to the 
      #   axis of the input. The number of elements in :real_real_kind must be equal to
      #   :dim. Can accept one of the inputs specified in REAL_REAL_FFT_KINDS_HASH.
      #   @see REAL_REAL_FFT_KINDS_HASH
      #   @see http://www.fftw.org/fftw3_doc/Real_002dto_002dReal-Transform-Kinds.html#Real_002dto_002dReal-Transform-Kinds
      # @example Create a plan for a basic 1D FFT and execute it.
      #   input = NMatrix.new([10],
      #     [
      #       Complex(9.32,0), Complex(44,0), Complex(125,0), Complex(34,0),
      #       Complex(31,0),   Complex(44,0), Complex(12,0),  Complex(1,0),
      #       Complex(53.23,0),Complex(-23.23,0),
      #     ], dtype: :complex128)
      #   plan = NMatrix::FFTW::Plan.new(10)
      #   plan.set_input input
      #   plan.execute
      #   print plan.output
      def initialize shape, opts={}
        verify_opts opts
        opts = {
          dim: 1,
          flags: :estimate,
          direction: :forward,
          type: :complex_complex
        }.merge(opts)

        @type      = opts[:type]
        @dim       = opts[:dim]
        @direction = opts[:direction]
        @shape     = shape.is_a?(Array) ? shape : [shape]
        @size      = @shape[0...@dim].inject(:*)
        @flags     = opts[:flags].is_a?(Array) ? opts[:flags] : [opts[:flags]]
        @real_real_kind    = opts[:real_real_kind]

        raise ArgumentError, ":real_real_kind option must be specified for :real_real type transforms" if
          @real_real_kind.nil? and @type == :real_real

        raise ArgumentError, "Specify kind of transform of each axis of input." if
          @real_real_kind and @real_real_kind.size != @dim

        raise ArgumentError, "dim (#{@dim}) cannot be more than size of shape #{@shape.size}" if
          @dim > @shape.size

        @plan_data = c_create_plan(@shape, @size, @dim, 
          combine_flags(@flags), FFT_DIRECTION_HASH[@direction], 
          DATA_TYPE_HASH[@type], encoded_rr_kind)
      end

      # Set input for the planned DFT.
      # @param [NMatrix] ip An NMatrix specifying the input to the FFT routine.
      #   The data type of the NMatrix must be either :complex128 or :float64
      #   depending on the type of FFT that has been planned. Size must be same
      #   as the size of the planned routine.
      # @raise [ArgumentError] if the input has any storage apart from :dense
      #   or if size/data type of the planned transform and the input matrix
      #   don't match.
      def set_input ip
        raise ArgumentError, "stype must be dense." if ip.stype != :dense
        raise ArgumentError, "size of input (#{ip.size}) cannot be greater than planned input size #{@size}" if
          ip.size != @size
        
        case @type
        when :complex_complex, :complex_real
          raise ArgumentError, "dtype must be complex128." if ip.dtype != :complex128
        when :real_complex, :real_real
          raise ArgumentError, "dtype must be float64." if ip.dtype != :float64
        else
          raise "Invalid type #{@type}"
        end

        @input = ip
        c_set_input(ip, @plan_data, DATA_TYPE_HASH[@type])
      end

      # Execute the DFT with the set plan.
      # @return [TrueClass] If all goes well and the fourier transform has been
      #   sucessfully computed, 'true' will be returned and you can access the
      #   computed output from the NMatrix::FFTW::Plan#output accessor.
      def execute
        @output = 
        case @type
        when :complex_complex
          @input.clone_structure        
        when :real_complex
          NMatrix.new([@input.size/2 + 1], dtype: :complex128)
        when :complex_real, :real_real
          NMatrix.new([@input.size], dtype: :float64)
        else
          raise TypeError, "Invalid type #{@type}"
        end

        c_execute(@output, @plan_data, DATA_TYPE_HASH[@type])
      end
     private

      # Combine flags received from the user (Symbols) into their respective
      # numeric equivalents and then 'OR' (|) all of them so the resulting number
      # can be passed directly to the FFTW planner function.
      def combine_flags flgs
        temp = 0
        flgs.each do |f|
          temp |= FLAG_VALUE_HASH[f]
        end
        temp
      end

      # Verify options passed into the constructor to make sure that no invalid
      # options have been passed.
      def verify_opts opts
        unless (opts.keys - VALID_OPTS).empty?
          raise ArgumentError, "#{opts.keys - VALID_OPTS} are invalid opts."
        end
      end

      # Get the numerical equivalents of the kind of real-real FFT to be computed.
      def encoded_rr_kind
        return @real_real_kind.map { |e| REAL_REAL_FFT_KINDS_HASH[e] } if @real_real_kind
      end
    end
  end
end