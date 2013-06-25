# -*- ruby -*-

require 'rubygems'
require 'bundler'
begin
  Bundler.setup(:default, :development)
rescue Bundler::BundlerError => e
  $stderr.puts e.message
  $stderr.puts "Run `bundle install` to install missing gems"
  exit e.status_code
end

require 'rake'
require "rake/extensiontask"
Rake::ExtensionTask.new do |ext|
    ext.name = 'nmatrix'          
    ext.ext_dir = 'ext/nmatrix' 
    ext.lib_dir = 'lib/'             
    ext.source_pattern = "**/*.{c,cpp, h}" 
end

gemspec = eval(IO.read("nmatrix.gemspec"))

require "rake/gempackagetask"
Rake::GemPackageTask.new(gemspec).define

desc "install the gem locally"
task :install => [:package] do
  sh %{gem install pkg/nmatrix-#{NMatrix::VERSION}.gem}
end

require 'rspec/core/rake_task'
require 'rspec/core'
require 'rspec/core/rake_task'
RSpec::Core::RakeTask.new(:spec) do |spec|
  spec.pattern = FileList['spec/**/*_spec.rb']
end

BASEDIR = Pathname( __FILE__ ).dirname.relative_path_from( Pathname.pwd )
SPECDIR = BASEDIR + 'spec'

VALGRIND_OPTIONS = [
        "--tool=memcheck",
        "--leak-check=yes",
        "--num-callers=15",
        #"--error-limit=no",
        "--partial-loads-ok=yes",
        "--undef-value-errors=no" #,
        #"--dsymutil=yes"
]
VALGRIND_MEMORYFILL_OPTIONS = [
        "--freelist-vol=100000000",
        "--malloc-fill=6D",
        "--free-fill=66 ",
]

GDB_OPTIONS = []


RSpec::Core::RakeTask.new(:spec)

task :console do |task|
  cmd = [ 'irb', "-r './lib/nmatrix.rb'" ]
  run *cmd
end

task :pry do |task|
  cmd = [ 'pry', "-r './lib/nmatrix.rb'" ]
  run *cmd
end

#namespace :console do
#  CONSOLE_CMD = ['irb', "-r './lib/nmatrix.rb'"]
#  desc "Run console under GDB."
#  task :gdb => [ :compile ] do |task|
#          cmd = [ 'gdb' ] + GDB_OPTIONS
#          cmd += [ '--args' ]
#          cmd += CONSOLE_CMD
#          run( *cmd )
#  end
#
#  desc "Run console under Valgrind."
#  task :valgrind => [ :compile ] do |task|
#          cmd = [ 'valgrind' ] + VALGRIND_OPTIONS
#          cmd += CONSOLE_CMD
#          run( *cmd )
#  end
#end

task :default => :spec

def run *cmd
  sh(cmd.join(" "))
end

namespace :spec do
  # partial-loads-ok and undef-value-errors necessary to ignore
  # spurious (and eminently ignorable) warnings from the ruby
  # interpreter

  RSPEC_CMD = [ 'ruby', '-S', 'rspec', '-Ilib:ext', SPECDIR ]

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
end


LEAKCHECK_CMD = [ 'ruby', '-Ilib:ext', "#{SPECDIR}/leakcheck.rb" ]

desc "Run leakcheck script."
task :leakcheck => [ :compile ] do |task|
  cmd = [ 'valgrind' ] + VALGRIND_OPTIONS
  cmd += LEAKCHECK_CMD
  run( *cmd )
end

namespace :clean do
  task :so do |task|
    tmp_path = "tmp/#{RUBY_PLATFORM}/nmatrix/#{RUBY_VERSION}"
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

# vim: syntax=ruby
