lib = File.expand_path('../lib/', __FILE__)
$:.unshift lib unless $:.include?(lib)

require 'nmatrix/version'

Gem::Specification.new do |gem|
  gem.name = "nmatrix"
  gem.version = NMatrix::VERSION::STRING
  gem.summary = "NMatrix is an experimental linear algebra library for Ruby, written mostly in C."
  gem.description = "NMatrix is an experimental linear algebra library for Ruby, written mostly in C."
  gem.homepage = 'http://sciruby.com'
  gem.authors = ['John Woods', 'Chris Wailes', 'Aleksey Timin']
  gem.email =  ['john.o.woods@gmail.com']
  gem.license = 'BSD 2-clause'
  gem.post_install_message = <<-EOF
***********************************************************
Welcome to SciRuby: Tools for Scientific Computing in Ruby!

NMatrix requires a C compiler, and has been tested only
with GCC 4.6+. We are happy to accept contributions
which improve the portability of this project.

Also required is ATLAS. Most Linux distributions and Mac
versions include ATLAS, but you may wish to compile it
yourself. The Ubuntu/Debian apt package for ATLAS WILL
NOT WORK with NMatrix if LAPACK is also installed.

More explicit instructions for NMatrix and SciRuby should
be available on the SciRuby website, sciruby.com, or
through our mailing list (which can be found on our web-
site).

Thanks for trying out NMatrix! Happy coding!

***********************************************************
EOF

  gem.files         = `git ls-files`.split("\n")
  gem.test_files    = `git ls-files -- {test,spec,features}/*`.split("\n")
  gem.executables   = `git ls-files -- bin/*`.split("\n").map{ |f| File.basename(f) }
  gem.extensions = ['ext/nmatrix/extconf.rb']
  gem.require_paths = ["lib"]

  gem.required_ruby_version = '>= 1.9'

  gem.add_dependency 'rdoc', '~>4.0', '>=4.0.1'
  gem.add_dependency 'packable', '>=1.3.5'
  gem.add_development_dependency 'rake', '~>10.3'
  gem.add_development_dependency 'bundler', '~>1.6'
  gem.add_development_dependency 'rspec', '~>2.14'
  gem.add_development_dependency 'rspec-longrun', '~>1.0'
  gem.add_development_dependency 'pry', '~>0.10'
  gem.add_development_dependency 'rake-compiler', '~>0.8'
end

