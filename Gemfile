source "https://rubygems.org"

# main gemspec
gemspec name: "nmatrix"

# plugin gemspecs
Dir["nmatrix-*.gemspec"].each do |gemspec_file|
  plugin_name = gemspec_file.match(/(nmatrix-.*)\.gemspec/)[1]
  gemspec(name: plugin_name, development_group: :plugin)
end
