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
# == homogeneous.rb
#
# This file contains optional shortcuts for generating homogeneous
# transformations.
#
#++

class NMatrix
  class << self
    #
    # call-seq:
    #     x_rotation(angle_in_radians) -> NMatrix
    #     x_rotation(angle_in_radians, dtype: dtype) -> NMatrix
    #     y_rotation(angle_in_radians) -> NMatrix
    #     y_rotation(angle_in_radians, dtype: dtype) -> NMatrix
    #     z_rotation(angle_in_radians) -> NMatrix
    #     z_rotation(angle_in_radians, dtype: dtype) -> NMatrix
    #
    # Generate a 4x4 homogeneous transformation matrix representing a rotation
    # about the x, y, or z axis respectively.
    #
    # * *Arguments* :
    #   - +angle_in_radians+ -> The angle of rotation in radians.
    #   - +dtype+ -> (optional) Default is +:float64+
    # * *Returns* :
    #   - A homogeneous transformation matrix consisting of a single rotation.
    #
    # Examples:
    #
    #    NMatrix.x_rotation(Math::PI.quo(6)) # =>
    #                                              1.0      0.0       0.0       0.0
    #                                              0.0      0.866025 -0.499999  0.0
    #                                              0.0      0.499999  0.866025  0.0
    #                                              0.0      0.0       0.0       1.0
    #
    #
    #    NMatrix.x_rotation(Math::PI.quo(6), dtype: :float32) # =>
    #                                              1.0      0.0       0.0       0.0
    #                                              0.0      0.866025 -0.5       0.0
    #                                              0.0      0.5       0.866025  0.0
    #                                              0.0      0.0       0.0       1.0
    #
    def x_rotation angle_in_radians, opts={}
      c = Math.cos(angle_in_radians)
      s = Math.sin(angle_in_radians)
      NMatrix.new(4, [1.0, 0.0, 0.0, 0.0,
                      0.0, c,   -s,  0.0,
                      0.0, s,    c,  0.0,
                      0.0, 0.0, 0.0, 1.0], {dtype: :float64}.merge(opts))
    end

    def y_rotation angle_in_radians, opts={}
      c = Math.cos(angle_in_radians)
      s = Math.sin(angle_in_radians)
      NMatrix.new(4, [ c,  0.0,  s,  0.0,
                      0.0, 1.0, 0.0, 0.0,
                      -s,  0.0,  c,  0.0,
                      0.0, 0.0, 0.0, 1.0], {dtype: :float64}.merge(opts))
    end

    def z_rotation angle_in_radians, opts={}
      c = Math.cos(angle_in_radians)
      s = Math.sin(angle_in_radians)
      NMatrix.new(4, [ c,  -s,  0.0, 0.0,
                       s,   c,  0.0, 0.0,
                      0.0, 0.0, 1.0, 0.0,
                      0.0, 0.0, 0.0, 1.0], {dtype: :float64}.merge(opts))
    end


    #
    # call-seq:
    #     translation(x, y, z) -> NMatrix
    #     translation([x,y,z]) -> NMatrix
    #     translation(translation_matrix) -> NMatrix
    #     translation(translation_matrix) -> NMatrix
    #     translation(translation, dtype: dtype) -> NMatrix
    #     translation(x, y, z, dtype: dtype) -> NMatrix
    #
    # Generate a 4x4 homogeneous transformation matrix representing a translation.
    #
    # * *Returns* :
    #   - A homogeneous transformation matrix consisting of a translation.
    #
    # Examples:
    #
    #    NMatrix.translation(4.0,5.0,6.0) # =>
    #                                          1.0   0.0   0.0   4.0
    #                                          0.0   1.0   0.0   5.0
    #                                          0.0   0.0   1.0   6.0
    #                                          0.0   0.0   0.0   1.0
    #
    #    NMatrix.translation(4.0,5.0,6.0, dtype: :int64) # =>
    #                                                         1  0  0  4
    #                                                         0  1  0  5
    #                                                         0  0  1  6
    #                                                         0  0  0  1
    #    NMatrix.translation(4,5,6) # =>
    #                                     1  0  0  4
    #                                     0  1  0  5
    #                                     0  0  1  6
    #                                     0  0  0  1
    #
    def translation *args
      xyz = args.shift if args.first.is_a?(NMatrix) || args.first.is_a?(Array)
      default_dtype = xyz.respond_to?(:dtype) ? xyz.dtype : NMatrix.guess_dtype(xyz)
      opts = {dtype: default_dtype}
      opts = opts.merge(args.pop) if args.size > 0 && args.last.is_a?(Hash)
      xyz ||= args

      n = if args.size > 0
        NMatrix.eye(4, opts)
      else
        NMatrix.eye(4, opts)
      end
      n[0..2,3] = xyz
      n
    end
  end
end