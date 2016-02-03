lib = File.expand_path('../lib/', __FILE__)
$:.unshift lib unless $:.include?(lib)

require 'nmatrix/version'

Gem::Specification.new do |gem|
  gem.name = "nmatrix-csv"
  gem.version = NMatrix::VERSION::STRING
  gem.summary = "C++ CSV reader for nmatrix"
  gem.description = "A csv reader implemented in C++ for nmatrix"
  gem.homepage = 'http://sciruby.com'
  gem.authors = ['Zhuanhao Wu']
  gem.email =  ['allen_ng@foxmail.com']
  gem.license = 'BSD 3-clause'

  gem.files         = ["lib/nmatrix/csv.rb"]
  gem.files         += `git ls-files -- ext/nmatrix_csv`.split("\n")
  gem.files         += `git ls-files -- ext/nmatrix | grep ".h$"`.split("\n") #need nmatrix header files to compile
  gem.test_files    = `git ls-files -- spec`.split("\n")
  gem.test_files    -= `git ls-files -- spec/plugins`.split("\n")
  gem.test_files    += `git ls-files -- spec/plugins/csv`.split("\n")
  gem.extensions = ['ext/nmatrix_csv/extconf.rb']
  gem.require_paths = ["lib"]

  gem.required_ruby_version = '>= 1.9'

  gem.add_dependency 'nmatrix', NMatrix::VERSION::STRING
end

