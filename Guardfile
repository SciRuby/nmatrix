guard 'rspec' do
  watch(%r{^spec/.+_spec\.rb$})
  watch(%r{^lib/(.+)\.rb$}) { "spec" }
  watch(%r{^ext/nmatrix/(.+)\.(c|h|rb)$}) { `rake compile | rake spec` }
  watch('spec/spec_helper.rb')  { "spec" }
end
