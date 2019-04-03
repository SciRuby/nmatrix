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
    def x_rotation angle_in_radians, opts = {}
      c = Math.cos(angle_in_radians)
      s = Math.sin(angle_in_radians)
      NMatrix.new(4, [1.0, 0.0, 0.0, 0.0,
                      0.0, c,   -s,  0.0,
                      0.0, s, c,  0.0,
                      0.0, 0.0, 0.0, 1.0,], {dtype: :float64}.merge(opts))
    end

    def y_rotation angle_in_radians, opts = {}
      c = Math.cos(angle_in_radians)
      s = Math.sin(angle_in_radians)
      NMatrix.new(4, [c, 0.0, s,  0.0,
                      0.0, 1.0, 0.0, 0.0,
                      -s,  0.0, c, 0.0,
                      0.0, 0.0, 0.0, 1.0,], {dtype: :float64}.merge(opts))
    end

    def z_rotation angle_in_radians, opts = {}
      c = Math.cos(angle_in_radians)
      s = Math.sin(angle_in_radians)
      NMatrix.new(4, [c, -s, 0.0, 0.0,
                      s, c,  0.0, 0.0,
                      0.0, 0.0, 1.0, 0.0,
                      0.0, 0.0, 0.0, 1.0,], {dtype: :float64}.merge(opts))
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
      n[0..2, 3] = xyz
      n
    end
  end

  #
  # call-seq:
  #     quaternion -> NMatrix
  #
  # Find the quaternion for a 3D rotation matrix.
  #
  # Code borrowed from: http://courses.cms.caltech.edu/cs171/quatut.pdf
  #
  # * *Returns* :
  #   - A length-4 NMatrix representing the corresponding quaternion.
  #
  # Examples:
  #
  #    n.quaternion # => [1, 0, 0, 0]
  #
  def quaternion
    raise(ShapeError, "Expected square matrix") if shape[0] != shape[1]
    raise(ShapeError, "Expected 3x3 rotation (or 4x4 homogeneous) matrix") if shape[0] > 4 || shape[0] < 3

    q = NMatrix.new([4], dtype: dtype == :float32 ? :float32 : :float64)
    rotation_trace = self[0, 0] + self[1, 1] + self[2, 2]
    if rotation_trace >= 0
      self_w = shape[0] == 4 ? self[3, 3] : 1.0
      root_of_homogeneous_trace = Math.sqrt(rotation_trace + self_w)
      q[0] = root_of_homogeneous_trace * 0.5
      s = 0.5 / root_of_homogeneous_trace
      q[1] = (self[2, 1] - self[1, 2]) * s
      q[2] = (self[0, 2] - self[2, 0]) * s
      q[3] = (self[1, 0] - self[0, 1]) * s
    else
      h = 0
      h = 1 if self[1, 1] > self[0, 0]
      h = 2 if self[2, 2] > self[h, h]

      case_macro = proc { |i, j, k, ii, jj, kk|
        qq = NMatrix.new([4], dtype: :float64)
        self_w = shape[0] == 4 ? self[3, 3] : 1.0
        s = Math.sqrt((self[ii, ii] - (self[jj, jj] + self[kk, kk])) + self_w)
        qq[i] = s * 0.5
        s = 0.5 / s
        qq[j] = (self[ii, jj] + self[jj, ii]) * s
        qq[k] = (self[kk, ii] + self[ii, kk]) * s
        qq[0] = (self[kk, jj] - self[jj, kk]) * s
        qq
      }

      case h
      when 0
        q = case_macro.call(1, 2, 3, 0, 1, 2)
      when 1
        q = case_macro.call(2, 3, 1, 1, 2, 0)
      when 2
        q = case_macro.call(3, 1, 2, 2, 0, 1)
      end

      self_w = shape[0] == 4 ? self[3, 3] : 1.0
      if self_w != 1
        s = 1.0 / Math.sqrt(self_w)
        q[0] *= s
        q[1] *= s
        q[2] *= s
        q[3] *= s
      end
    end

    q
  end

  #
  # call-seq:
  #     angle_vector -> [angle, about_vector]
  #
  # Find the angle vector for a quaternion. Assumes the quaternion has unit length.
  #
  # Source: http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/
  #
  # * *Returns* :
  #   - An angle (in radians) describing the rotation about the +about_vector+.
  #   - A length-3 NMatrix representing the corresponding quaternion.
  #
  # Examples:
  #
  #    q.angle_vector # => [1, 0, 0, 0]
  #
  def angle_vector
    raise(ShapeError, "Expected length-4 vector or matrix (quaternion)") if shape[0] != 4
    raise("Expected unit quaternion") if self[0] > 1

    xyz = NMatrix.new([3], dtype: dtype)

    angle = 2 * Math.acos(self[0])
    s = Math.sqrt(1.0 - self[0] * self[0])

    xyz[0..2] = self[1..3]
    xyz /= s if s >= 0.001 # avoid divide by zero
    [angle, xyz]
  end
end
