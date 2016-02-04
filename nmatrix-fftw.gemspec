lib = File.expand_path('../lib/', __FILE__)
$:.unshift lib unless $:.include?(lib)

require 'nmatrix/version'

Gem::Specification.new do |gem|
  gem.name = "nmatrix-fftw"
  gem.version = NMatrix::VERSION::STRING
  gem.summary = "FFTW backend for nmatrix"
  gem.description = "NMatrix extension for using fuctions provided by FFTW"
  gem.homepage = 'http://sciruby.com'
  gem.authors = ['Sameer Deshmukh', 'Magdalen Berns']
  gem.email =  ['sameer.deshmukh93@gmail.com', 'm.berns@thismagpie.com']
  gem.license = 'BSD 3-clause'

  gem.files         = ["lib/nmatrix/fftw.rb"]
  gem.files         += `git ls-files -- ext/nmatrix_fftw`.split("\n")
  gem.files         += `git ls-files -- ext/nmatrix | grep ".h$"`.split("\n") #need nmatrix header files to compile
  gem.test_files    = `git ls-files -- spec`.split("\n")
  gem.test_files    -= `git ls-files -- spec/plugins`.split("\n")
  gem.test_files    += `git ls-files -- spec/plugins/fftw`.split("\n")
  gem.extensions = ['ext/nmatrix_fftw/extconf.rb']
  gem.require_paths = ["lib"]

  gem.required_ruby_version = '>= 1.9'

  gem.add_dependency 'nmatrix', NMatrix::VERSION::STRING
end

