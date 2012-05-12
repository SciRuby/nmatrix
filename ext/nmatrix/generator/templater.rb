module Generator

  # Represents a single output C file, generated from a series of templates.
  #
  # Example usage and instructions
  # ==============================
  #
  #    Generator::Templater.new('blas.c', :in => 'dense', :boilerplate => 'blas_header') do |c|
  #      c.template 'rationalmath', :TYPE => Generator::RATIONAL_DTYPES
  #      c.template 'complexmath', :in => 'yale', :TYPE => Generator::COMPLEX_DTYPES
  #      c.template %w{gemm gemv}, :TYPE => Generator::NONBLAS_DTYPES
  #      c.update_header 'nmatrix.h'
  #    end
  #
  # In the example above, the templates used will be dense/rationalmath.template.c, dense/gemm.template.c,
  # dense/gemv.template.c, and yale/complexmath.template.c.
  #
  # Options which look like constants (e.g., :TYPE and :INT) specify the template types. The corresponding hash values
  # should be arrays of DTypeInfo objects.
  #
  # Any options taken by template can go in the call to Templater.new, and then need not be specified for individual
  # templates.
  #
  # You may also do:
  #
  #     c.boilerplate 'blas_header', :in => 'yale'
  #
  # as long as you do it before the template calls.
  #
  class Templater

    # Create a C file from a template. You must supply the name of the C file to write.
    #
    # If you wish to pass options to add by default, you may pass them in the second argument.
    # You may also supply a :boilerplate option, giving the filename of the output's boilerplate, if you so desire.
    def initialize output_filename, default_options = {}
      @output_filename = output_filename
      @boiler_filename = default_options.delete(:boilerplate)
      @header_filename = default_options.delete(:header)
      @default_options = default_options
      @declarations    = []

      @begun = false
      unless @boiler_filename.nil?
        boilerplate(@boiler_filename, default_options)
        @begun = true
      end

      # This yield assigns @declarations.
      yield self

      update_header(@update_header) unless @update_header.nil?
    end

    def boilerplate name, options = {}
      raise("please provide boilerplate before templates") if @begun
      add_internal(name, false, options)
    end

    # Add a .template.c file to the current output. Takes templates as mandatory argument (e.g., for numbmm.template.c
    # and transp.template.c, give it %w{numbmm transp}).
    #
    # Options
    # =======
    #  * :in  :: subdirectory to search in
    def template names, options = {}
      names = [names] if names.is_a?(String)
      @declarations = @declarations.concat(add_internal(names, true, options))
    end

    # Modify some .h file by adding the templated function declarations to it. This function is called automatically
    # if you set :header in the creation of the Templater, but you may also call it manually -- as long as you do so
    # after the processing of the templates whose prototypes you want to include.
    #
    # This takes no options, and ignores default_options. It expects the header to be in the root source directory.
    def update_header name
      from_filename = "#{name}.inc.h"
      filename      = "#{name}.h"

      # Find the current version of the inc file. It may have already been written in the build directory, in which case
      # we should start with that one rather than the true original.
      orig =
      if File.exists?(filename)
        `mv #{filename} #{from_filename}`
        File.new(from_filename, "r") # read the temporary version in the src directory
      else
        File.new(relative_ext_path(from_filename), "r") # read the one in the src directory
      end

      out = File.new(filename, "w") # write the one in the build directory

      while line = orig.gets
        line.chomp!
        if line.include?("%%INSERT_TEMPLATER_DECLARATIONS%%")
          Generator.process_declarations(@declarations).each { |prototype| out.puts prototype }
        end

        # Output the cue even if we added declarations, because later instantiations may also add declarations.
        out.puts line
      end
    end

  protected

    def add_internal names, as_template, options
      opts = options.dup
      dir = opts.delete(:in) || @default_options[:in]

      if as_template
        template_hash = @default_options.merge(opts).select { |k,v| k.to_s =~ /^[A-Z_0-9]+$/}
        apply_each_template_type(@output_filename, names, dir, template_hash)
      else
        # Pretend names is just one name for this part:
        template_path = File.join(relative_ext_path(dir), "#{names}.template.c")
        `cat #{template_path} > ./#{@output_filename}`
        []
      end
    end

    def relative_ext_path subdir=nil
      subdir.nil? || subdir.size == 0 ? File.join($RELATIVE_PATH, "ext", "nmatrix") : File.join($RELATIVE_PATH, "ext", "nmatrix", subdir)
    end

    # Returns an array of C function declarations, which can then be put in nmatrix.h
    def apply_each_template_type output_path, names, dir, type_hash
      raise(NotImplementedError, "Too many keys, please fix") if type_hash.keys.size > 2

      template_path = !dir.nil? && dir.size > 0 ? File.join(relative_ext_path, dir) : File.join($RELATIVE_PATH, "ext", "nmatrix")

      declarations = []

      a = type_hash.keys[0]
      b = type_hash.keys[1] if type_hash.size > 1

      raise(ArgumentError, "no template keys found") if !defined?(a) || a.nil?

      type_hash[a].each do |type_a|

        # Currently allows just one or two types. TODO: Fix this so it allows as many types as we could want. Low priority.
        if defined?(b) && !b.nil?
          type_hash[b].each do |type_b|
            names.each do |name|
              full_template_path = File.join(template_path, "#{name}.template.c")
              raise(ArgumentError, "template not found at #{full_template_path} relative to #{Dir.pwd}") unless File.exists?(full_template_path)

              declarations = declarations.concat(Generator.template(full_template_path, output_path, {a => type_a, b => type_b}))
            end
          end
        else
          names.each do |name|

            full_template_path = File.join(template_path, "#{name}.template.c")
            raise(ArgumentError, "template not found at #{full_template_path} relative to #{Dir.pwd}") unless File.exists?(full_template_path)

            declarations = declarations.concat(Generator.template(full_template_path, output_path, {a => type_a}))
          end
        end

      end

      declarations
    end
  end
end
