lib = File.expand_path("../lib/", __FILE__)
$:.unshift lib unless $:.include?(lib)

require "nmatrix/version"

Gem::Specification.new do |gem|
  gem.name = "nmatrix-lapacke"
  gem.version = NMatrix::VERSION::STRING
  gem.summary = "general LAPACK backend for nmatrix using LAPACKE interface"
  gem.description = "For using linear algebra fuctions provided by LAPACK and BLAS"
  gem.homepage = "http://sciruby.com"
  gem.authors = ["Will Levine", "John Woods"]
  gem.email = ["john.o.woods@gmail.com"]
  gem.license = "BSD-3-Clause"

  gem.files = ["lib/nmatrix/lapacke.rb", "lib/nmatrix/lapack_ext_common.rb"]
  gem.files         += `git ls-files -- ext/nmatrix_lapacke`.split("\n")
  gem.files         += `git ls-files -- ext/nmatrix | grep ".h$"`.split("\n") # need nmatrix header files to compile
  gem.test_files = `git ls-files -- spec`.split("\n")
  gem.test_files    -= `git ls-files -- spec/plugins`.split("\n")
  gem.test_files    += `git ls-files -- spec/plugins/lapacke`.split("\n")
  gem.extensions = ["ext/nmatrix_lapacke/extconf.rb"]
  gem.require_paths = ["lib"]

  gem.required_ruby_version = ">= 1.9"

  gem.add_dependency "nmatrix", NMatrix::VERSION::STRING
end
