# -*- ruby -*-

require 'rubygems'
require 'rubygems/package_task'
require 'bundler'

#Specify plugins to build on the command line like:
#rake whatever nmatrix_plugins=atlas,lapacke
#or
#rake whatever nmatrix_plugins=all
#If you want to build *only* plugins and not the core nmatrix gem:
#rake whatever nmatrix_plugins=all nmatrix_core=false
if ENV["nmatrix_plugins"] == "all"
  gemspecs = Dir["*.gemspec"]
else
  plugins = []
  plugins = ENV["nmatrix_plugins"].split(",") if ENV["nmatrix_plugins"]
  gemspecs = ["nmatrix.gemspec"] #always include the main nmatrix gem
  plugins.each do |plugin|
    gemspecs << "nmatrix-#{plugin}.gemspec"
  end
end
if ENV["nmatrix_core"] == "false"
  gemspecs -= ["nmatrix.gemspec"]
end
gemspecs.map! { |gemspec| eval(IO.read(gemspec)) }

begin
  Bundler.setup(:default, :development)
rescue Bundler::BundlerError => e
  $stderr.puts e.message
  $stderr.puts "Run `bundle install` to install missing gems"
  exit e.status_code
end

desc "Build and install into system gems."
task :install => :repackage do
  gemspecs.each do |gemspec|
    gem_file = "pkg/#{gemspec.name}-#{gemspec.version}.gem"
    system "gem install '#{gem_file}'"
  end
end

require 'rake'
require "rake/extensiontask"

gemspecs.each do |gemspec|
  next unless gemspec.extensions
  gemspec.extensions.each do |extconf|
    ext_name = extconf.match(/ext\/(.*)\/extconf\.rb/)[1]
    Rake::ExtensionTask.new do |ext|
      ext.name = ext_name
      ext.ext_dir = "ext/#{ext_name}"
      ext.lib_dir = 'lib/'
      ext.source_pattern = "**/*.{c,cpp,h}"
    end
  end
end

gemspecs.each do |gemspec|
  Gem::PackageTask.new(gemspec).define
end

require 'rspec/core/rake_task'
require 'rspec/core'
namespace :spec do
  #We need a separate rake task for each plugin, rather than one big task that
  #runs all of the specs. This is because there's no way to tell rspec
  #to run the specs in a certain order with (say) "nmatrix/atlas" require'd
  #for some of the specs, but not for others, without splitting them up like
  #this.
  spec_tasks = []
  gemspecs.each do |gemspec|
    test_files = gemspec.test_files
    test_files.keep_if { |file| file =~ /_spec\.rb$/ }
    next if test_files.empty?
    spec_tasks << gemspec.name
    RSpec::Core::RakeTask.new(gemspec.name) do |spec|
      spec.pattern = FileList.new(test_files)
    end
  end
  task :all => spec_tasks
end



task :spec => "spec:all"

BASEDIR = Pathname( __FILE__ ).dirname.relative_path_from( Pathname.pwd )
SPECDIR = BASEDIR + 'spec'

VALGRIND_OPTIONS = [
    "--tool=memcheck",
    #"--leak-check=yes",
    "--num-callers=15",
    #"--error-limit=no",
    "--partial-loads-ok=yes",
    "--undef-value-errors=no" #,
    #"--dsymutil=yes"
]

CALLGRIND_OPTIONS = [
    "--tool=callgrind",
    "--dump-instr=yes",
    "--simulate-cache=yes",
    "--collect-jumps=yes"
]

VALGRIND_MEMORYFILL_OPTIONS = [
    "--freelist-vol=100000000",
    "--malloc-fill=6D",
    "--free-fill=66 ",
]

GDB_OPTIONS = []


task :console do |task|
  cmd = [ 'irb', "-r './lib/nmatrix.rb'" ]
  run *cmd
end

task :pry do |task|
  cmd = [ 'pry', "-r './lib/nmatrix.rb'" ]
  run *cmd
end

namespace :pry do
  task :valgrind => [ :compile ] do |task|
    cmd  = [ 'valgrind' ] + VALGRIND_OPTIONS
    cmd += ['ruby', '-Ilib:ext', "-r './lib/nmatrix.rb'", "-r 'pry'", "-e 'binding.pry'"]
    run *cmd
  end
end

namespace :console do
  CONSOLE_CMD = ['irb', "-r './lib/nmatrix.rb'"]
  desc "Run console under GDB."
  task :gdb => [ :compile ] do |task|
          cmd = [ 'gdb' ] + GDB_OPTIONS
          cmd += [ '--args' ]
          cmd += CONSOLE_CMD
          run( *cmd )
  end

  desc "Run console under Valgrind."
  task :valgrind => [ :compile ] do |task|
          cmd = [ 'valgrind' ] + VALGRIND_OPTIONS
          cmd += CONSOLE_CMD
          run( *cmd )
  end
end

task :default => :spec

def run *cmd
  sh(cmd.join(" "))
end

namespace :spec do
  # partial-loads-ok and undef-value-errors necessary to ignore
  # spurious (and eminently ignorable) warnings from the ruby
  # interpreter

  RSPEC_CMD = [ 'ruby', '-S', 'rspec', '-Ilib:ext', SPECDIR.to_s ]

  #desc "Run the spec for generator.rb"
  #task :generator do |task|
  #  run 'rspec spec/generator_spec.rb'
  #end

  desc "Run specs under GDB."
  task :gdb => [ :compile ] do |task|
          cmd = [ 'gdb' ] + GDB_OPTIONS
    cmd += [ '--args' ]
    cmd += RSPEC_CMD
    run( *cmd )
  end

  desc "Run specs under cgdb."
  task :cgdb => [ :compile ] do |task|
    cmd = [ 'cgdb' ] + GDB_OPTIONS
    cmd += [ '--args' ]
    cmd += RSPEC_CMD
    run( *cmd )
  end

  desc "Run specs under Valgrind."
  task :valgrind => [ :compile ] do |task|
    cmd = [ 'valgrind' ] + VALGRIND_OPTIONS
    cmd += RSPEC_CMD
    run( *cmd )
  end

  desc "Run specs under Callgrind."
  task :callgrind => [ :compile ] do |task|
    cmd = [ 'valgrind' ] + CALLGRIND_OPTIONS
    cmd += RSPEC_CMD
    run( *cmd )
  end

end


LEAKCHECK_CMD = [ 'ruby', '-Ilib:ext', "#{SPECDIR}/leakcheck.rb" ]


desc "Run leakcheck script."
task :leakcheck => [ :compile ] do |task|
  cmd = [ 'valgrind' ] + VALGRIND_OPTIONS
  cmd += LEAKCHECK_CMD
  run( *cmd )
end

namespace :clean do
  #the generated Makefile doesn't have a soclean target, should this be removed?
  task :so do |task|
    gemspecs.each do |gemspec|
      next unless gemspec.extensions
      gemspec.extensions.each do |extconf|
        ext_name = extconf.match(/ext\/(.*)\/extconf\.rb/)[1]
        tmp_path = "tmp/#{RUBY_PLATFORM}/#{ext_name}/#{RUBY_VERSION}"
        chdir tmp_path do
          if RUBY_PLATFORM =~ /mswin/
            `nmake soclean`
          else
            mkcmd = ENV['MAKE'] || %w[gmake make].find { |c| system("#{c} -v >> /dev/null 2>&1") }
            `#{mkcmd} soclean`
          end
        end
      end
    end
  end
end


desc "Check the manifest for correctness"
task :check_manifest do |task|
  manifest_files  = File.read("Manifest.txt").split

  git_files       = `git ls-files |grep -v 'spec/'`.split
  ignore_files    = %w{.gitignore .rspec ext/nmatrix/binary_format.txt ext/nmatrix/ttable_helper.rb scripts/mac-brew-gcc.sh}

  possible_files  = git_files - ignore_files

  missing_files   = possible_files - manifest_files
  extra_files     = manifest_files - possible_files

  unless missing_files.empty?
    STDERR.puts "The following files are in the git repo but not the Manifest:"
    missing_files.each { |f| STDERR.puts " -- #{f}"}
  end

  unless extra_files.empty?
    STDERR.puts "The following files are in the Manifest but may not be necessary:"
    extra_files.each { |f| STDERR.puts " -- #{f}"}
  end

  if extra_files.empty? && missing_files.empty?
    STDERR.puts "Manifest looks good!"
  end

end

require "rdoc/task"
#separate out docs for plugins?
RDoc::Task.new do |rdoc|
  rdoc.main = "README.rdoc"
  rdoc.rdoc_files.include(%w{README.rdoc History.txt LICENSE.txt CONTRIBUTING.md lib ext})
  rdoc.options << "--exclude=ext/nmatrix/extconf.rb"
  rdoc.options << "--exclude=ext/nmatrix_atlas/extconf.rb"
  rdoc.options << "--exclude=ext/nmatrix/ttable_helper.rb"
  rdoc.options << "--exclude=lib/nmatrix/rspec.rb"
end

# jruby tasks

namespace :jruby do

  PROJECT_DIR = File.expand_path(".",Dir.pwd)

  BUILD_DIR = "build"
  CLASSES_DIR = "../build/classes"
  TEST_CLASSES_DIR = "build/testClasses"

  JRUBY_DIR = "#{PROJECT_DIR}/ext/nmatrix_java"
  VENDOR_DIR = "#{JRUBY_DIR}/vendor"
  TARGET_DIR = "#{JRUBY_DIR}/target"

  jars = Dir["#{VENDOR_DIR}/*.jar"]

  desc 'Compile java classes'
  task :javac do
    unless RUBY_PLATFORM == 'java'
      abort 'Please run with JRuby'
    end
    sh "mkdir -p #{JRUBY_DIR}/build/classes"
    Dir.chdir("#{JRUBY_DIR}/nmatrix")
    classes    = Dir['**/*.java']
    sh "javac -classpath #{jars.join(':')} -d #{CLASSES_DIR} #{classes.join(' ')}"
  end

  desc 'Package java classes in a jar file'
  task :jar do
    unless RUBY_PLATFORM == 'java'
      abort 'Please run with JRuby'
    end
    sh "mkdir -p #{TARGET_DIR}"
    Dir.chdir("#{JRUBY_DIR}/build/classes")
    classes = Dir['**/*.class']
    sh "jar -cf #{TARGET_DIR}/nmatrix.jar #{classes.join(' ')}"
  end

  task :all => [:javac, :jar]
end

desc "Compile java classes and Package them in a jar file"
task :jruby => 'jruby:all'

namespace :travis do
  task :env do
    puts "\n# Build environment:"
    %w[
      CC CXX
      USE_ATLAS USE_OPENBLAS USE_REF NO_EXTERNAL_LIB
      TRAVIS_OS_NAME TRAVIS_BRANCH TRAVIS_COMMIT TRAVIS_PULL_REQUEST
    ].each do |name|
      puts "- #{name}: #{ENV[name]}"
    end

    require 'rbconfig'
    puts "\n# RbConfig::MAKEFILE_CONFIG values:"
    %w[
      CC CXX CPPFLAGS CFLAGS CXXFLAGS
    ].each do |name|
      puts "- #{name}: #{RbConfig::MAKEFILE_CONFIG[name]}"
    end

    cc = RbConfig::MAKEFILE_CONFIG['CC']
    puts "\n$ #{cc} -v\n#{`#{cc} -v 2>&1`}"
  end
end

# vim: syntax=ruby
