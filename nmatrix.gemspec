lib = File.expand_path('../lib/', __FILE__)
$:.unshift lib unless $:.include?(lib)

require 'nmatrix/version'

#get files that are used by plugins rather than the main nmatrix gem
plugin_files = []
Dir["nmatrix-*.gemspec"].each do |gemspec_file|
  gemspec = eval(File.read(gemspec_file))
  plugin_files += gemspec.files
end
plugin_lib_files = plugin_files.select { |file| file.match(/^lib\//) }

Gem::Specification.new do |gem|
  gem.name = "nmatrix"
  gem.version = NMatrix::VERSION::STRING
  gem.summary = "NMatrix is a linear algebra library for Ruby"
  gem.description = "NMatrix is a linear algebra library for Ruby, written mostly in C and C++."
  gem.homepage = 'http://sciruby.com'
  gem.authors = ['John Woods', 'Chris Wailes', 'Aleksey Timin']
  gem.email =  ['john.o.woods@gmail.com']
  gem.license = 'BSD-3-Clause'
  gem.post_install_message = <<-EOF
***********************************************************
Welcome to SciRuby: Tools for Scientific Computing in Ruby!

NMatrix requires a C/C++ compiler. Clang and GCC are 
recommended. JRuby support is experimental, and requires
Java.

If you are upgrading from NMatrix 0.1.0 and rely on
ATLAS features, please check the README.

Faster matrix calculations and more advanced linear
algebra features are available by installing either
the nmatrix-atlas or nmatrix-lapacke plugins.

More explicit instructions for NMatrix and SciRuby should
be available on the SciRuby website, sciruby.com, or
through our mailing list (which can be found on our web-
site).

Thanks for trying out NMatrix! Happy coding!

***********************************************************
EOF

  gem.files         = `git ls-files -- ext/nmatrix`.split("\n")
  gem.files         += `git ls-files -- lib`.split("\n")
  gem.files         -= plugin_lib_files
  gem.test_files    = `git ls-files -- spec`.split("\n")
  gem.test_files    -= `git ls-files -- spec/plugins`.split("\n")
  gem.extensions = ['ext/nmatrix/extconf.rb']
  gem.require_paths = ["lib"]

  gem.required_ruby_version = '>= 1.9'

  gem.add_dependency 'packable', '~> 1.3', '>= 1.3.5'
  gem.add_development_dependency 'bundler', '>=1.6'
  gem.add_development_dependency 'pry', '~>0.10'
  gem.add_development_dependency 'rake', '~>10.3'
  gem.add_development_dependency 'rake-compiler', '~>0.8'
  gem.add_development_dependency 'rdoc', '~>4.0', '>=4.0.1'
  gem.add_development_dependency 'rspec', '~>2.14'
  gem.add_development_dependency 'rspec-longrun', '~>1.0'
end

