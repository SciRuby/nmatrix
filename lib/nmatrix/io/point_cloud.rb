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
# == io/point_cloud.rb
#
# Point Cloud Library (PCL) PCD file IO functions.
#
#++

# Reader for Point Cloud Data (PCD) file format.
#
# The documentation of this format can be found in:
#
# http://pointclouds.org/documentation/tutorials/pcd_file_format.php
#
# Note that this implementation does not take the width or height parameters
# into account.
module NMatrix::IO::PointCloud

  # For UINT, just add 1 to the index.
  INT_DTYPE_BY_SIZE   = [:int8, :int8, :int16, :int32, :int64, :int64] #:nodoc:
  FLOAT_DTYPE_BY_SIZE = {4 => :float32, 8 => :float64} #:nodoc:

  class << self
    # call-seq:
    #     load(filename) -> NMatrix
    #
    # * *Arguments* :
    #   - +filename+ -> String giving the name of the file to be loaded.
    #
    # Load a Point Cloud Library PCD file as a matrix.
    def load(filename)
      MetaReader.new(filename).matrix
    end
  end

  class MetaReader #:nodoc:
    ENTRIES = [:version,  :fields,           :size,  :type,
               :count,  :width,  :height,  :viewpoint,  :points,  :data]
    ASSIGNS = [:version=, :fields=,          :size=, :type=,
               :count=, :width=, :height=, :viewpoint=, :points=, :data=]
    CONVERT = [:to_s,     :downcase_to_sym,  :to_i,  :downcase_to_sym,
      :to_i,   :to_i,   :to_i,    :to_f,       :to_i,    :downcase_to_sym]

    DTYPE_CONVERT = {:byte => :to_i, :int8 => :to_i, :int16 => :to_i,
           :int32 => :to_i, :float32 => :to_f, :float64 => :to_f}

    # For UINT, just add 1 to the index.
    INT_DTYPE_BY_SIZE   = {1 => :int8,    2 => :int16,   4 => :int32,
       8 => :int64,  16 => :int64}
    FLOAT_DTYPE_BY_SIZE = {1 => :float32, 2 => :float32, 4 => :float32,
       8 => :float64,16 => :float64}

    class << self

      # Given a type and a number of bytes, figure out an appropriate dtype
      def dtype_by_type_and_size t, s
        if t == :f
          FLOAT_DTYPE_BY_SIZE[s]
        elsif t == :u
          return :byte if s == 1
          INT_DTYPE_BY_SIZE[s*2]
        else
          INT_DTYPE_BY_SIZE[s]
        end
      end
    end

    # call-seq:
    #     PointCloudReader::MetaReader.new(filename) -> MetaReader
    #
    # * *Arguments* :
    #   - +filename+ -> String giving the name of the file to be loaded.
    # * *Raises* :
    #   - +NotImplementedError+ -> only ASCII supported currently
    #   - +IOError+ -> premature end of file
    #
    # Open a file and read the metadata at the top; then read the PCD into an
    # NMatrix.
    #
    # In addition to the fields in the PCD file, there will be at least one
    # additional attribute, :matrix, storing the data.
    def initialize filename
      f = File.new(filename, "r")

      ENTRIES.each.with_index do |entry,i|
        read_entry(f, entry, ASSIGNS[i], CONVERT[i])
      end

      raise(NotImplementedError, "only ASCII supported currently") \
       unless self.data.first == :ascii

      @matrix = NMatrix.new(self.shape, dtype: self.dtype)

      # Do we want to use to_i or to_f?
      convert = DTYPE_CONVERT[self.dtype]

      i = 0
      while line = f.gets
        @matrix[i,:*] = line.chomp.split.map { |f| f.send(convert) }
        i += 1
      end

      raise(IOError, "premature end of file") if i < self.points[0]

    end

    attr_accessor *ENTRIES
    attr_reader :matrix

  protected
    # Read the current entry of the header.
    def read_entry f, entry, assign=nil, convert=nil
      assign ||= (entry.to_s + "=").to_sym

      while line = f.gets
        next if line =~ /^\s*#/ # ignore comment lines
        line = line.chomp.split(/\s*#/)[0] # ignore the comments after any data

        # Split, remove the entry name, and convert to the correct type.
        self.send(assign,
                  line.split.tap { |t| t.shift }.map do |f|
                    if convert.nil?
                      f
                    elsif convert == :downcase_to_sym
                      f.downcase.to_sym
                    else
                      f.send(convert)
                    end
                  end)

        # We don't really want to loop.
        break
      end

      self.send(entry)
    end


    # Determine the dtype for a matrix based on the types and
    #  sizes given in the PCD.
    #  Call this only after read_entry has been called.
    def dtype
      @dtype ||= begin
        dtypes = self.type.map.with_index do |t,k|
          MetaReader.dtype_by_type_and_size(t, size[k])
        end.sort.uniq

        # This could probably save one comparison at most, but we assume that
        # worst case isn't going to happen very often.
        while dtypes.size > 1
          d = NMatrix.upcast(dtypes[0], dtypes[1])
          dtypes.shift
          dtypes[0] = d
        end

        dtypes[0]
      end
    end

    # Determine the shape of the matrix.
    def shape
      @shape ||= [
          self.points[0],
          self.fields.size
      ]
    end
  end
end
