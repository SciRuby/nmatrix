NMatrix is part of SciRuby, a collaborative effort to bring scientific computation to Ruby. If you want to help, please
 do so!

This guide covers ways in which you can contribute to the development of SciRuby and, more specifically, NMatrix.

## How to help

There are various ways to help NMatrix: bug reports, coding and documentation. All of them are important.

First, you can help implement new features or bug fixes. To do that, visit our
 [roadmap](https://github.com/SciRuby/nmatrix/wiki/Roadmap) or our [issue tracker][2]. If you find something that you
 want to work on, post it in the issue or on our [mailing list][1].

You need to send tests together with your code. No exceptions. You can ask for our opinion, but we won't accept patches
 without good spec coverage.

We use RSpec for testing. If you aren't familiar with it, there's a good
 [guide to better specs with RSpec](http://betterspecs.org/) that shows a bit of the syntax and how to use it properly.
 However, the best resource is probably the specs that already exist -- so just read them.

And don't forget to write documentation (we use RDoc). It's necessary to allow others to know what's available in the
 library. There's a section on it later in this guide.

We only accept bug reports and pull requests in GitHub. You'll need to create a new (free) account if you don't have one
 already. To learn how to create a pull request, please see
 [this guide on collaborating](https://help.github.com/categories/63/articles).

If you have a question about how to use NMatrix or SciRuby in general or a feature/change in mind, please ask the
 [sciruby-dev mailing list][1].

Thanks!

## Coding

To start helping with the code, you need to have all the dependencies in place:

- ATLAS and LAPACK
- GCC 4.3+
- git
- Ruby 1.9+
- `bundler` gem

Now, you need to clone the git repository:

    git clone git://github.com/SciRuby/nmatrix.git
		cd nmatrix
		bundle install
		rake compile
		rake spec

This will install all dependencies, compile the extension and run the specs. 

As of now (12/31/2012), there should be 25 specs failing, in elementwise\_spec, lapack\_spec and math\_spec. If you see
more than 25 or from different specs, please report on the [mailing list][1] or on the [issue tracker][2].

If everything's fine until now, you can create a new branch to work on your feature:

    git branch new-feature
    git checkout new-feature

Before commiting any code, please read our
[Contributor Agreement](http://github.com/SciRuby/sciruby/wiki/Contributor-Agreement).

## Documentation

There are two ways in which NMatrix is being documented: guides and comments, which are converted with RDoc into the
documentation seen in [sciruby.com](http://sciruby.com).

If you want to write a guide on how to use NMatrix to solve some problem or simply showing how to use one of its
features, write it as a wiki page and send an e-mail on the [mailing list][1]. We're working to improve this process.

If you aren't familiar with RDoc syntax,
[this is the official documentation](http://docs.seattlerb.org/rdoc/RDoc/Markup.html).

## Conclusion

This guide was heavily based on the
[Contributing to Ruby on Rails guide](http://edgeguides.rubyonrails.org/contributing_to_ruby_on_rails.html).

[1]: https://groups.google.com/forum/?fromgroups#!forum/sciruby-dev
[2]: https://github.com/sciruby/nmatrix/issues?sort=created&state=open
