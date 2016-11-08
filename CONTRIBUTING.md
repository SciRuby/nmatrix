NMatrix is part of SciRuby, a collaborative effort to bring scientific
computation to Ruby. If you want to help, please do so!

This guide covers ways in which you can contribute to the development
of SciRuby and, more specifically, NMatrix.

## How to help

There are various ways to help NMatrix: bug reports, coding and
documentation. All of them are important.

First, you can help implement new features or bug fixes. To do that,
visit our [roadmap](https://github.com/SciRuby/nmatrix/wiki/Roadmap)
or our [issue tracker][2]. If you find something that you want to work
on, post it in the issue or on our [mailing list][1].

You need to send tests together with your code. No exceptions. You can
ask for our opinion, but we won't accept patches without good spec
coverage.

We use RSpec for testing. If you aren't familiar with it, there's a
good [guide to better specs with RSpec](http://betterspecs.org/) that
shows a bit of the syntax and how to use it properly.  However, the
best resource is probably the specs that already exist -- so just read
them.

And don't forget to write documentation (we use rdoc). It's necessary
to allow others to know what's available in the library. There's a
section on it later in this guide.

We only accept bug reports and pull requests in GitHub. You'll need to
create a new (free) account if you don't have one already. To learn
how to create a pull request, please see
[this guide on collaborating](https://help.github.com/categories/63/articles).

If you have a question about how to use NMatrix or SciRuby in general
or a feature/change in mind, please ask the
[sciruby-dev mailing list][1].

Thanks!

## Coding

To start helping with the code, you need to have all the dependencies in place:

- GCC 4.3+
- git
- Ruby 1.9+
- `bundler` gem
- ATLAS/LAPACKE/FFTW dependending on the plugin you want to change.

Now, you need to clone the git repository:

```bash
$ git clone git://github.com/SciRuby/nmatrix.git
$ cd nmatrix
$ bundle install
$ rake compile
$ rake spec
```

This will install all dependencies, compile the extension and run the
specs.

For **JRuby**

```bash
$ mkdir ext/nmatrix_java/vendor
Download commons_math.3.6.1 jar and place it in ext/nmatrix_java/vendor directory
$ mkdir -p ext/nmatrix_java/build/class
$ mkdir ext/nmatrix_java/target
$ rake jruby
```

If everything's fine until now, you can create a new branch to work on
your feature:

```bash
$ git branch new-feature
$ git checkout new-feature
```

Before commiting any code, please read our
[Contributor Agreement](http://github.com/SciRuby/sciruby/wiki/Contributor-Agreement).

### Guidelines for interfacing with C/C++

NMatrix uses a lot of C/C++ to speed up execution of processes and
give more control over data types, storage types, etc. Since we are
interfacing between two very different languages, things can get out
of hand pretty fast.

Please go thorough this before you create any C accessors:

* Perform all pre-computation error checking in Ruby.
* Curate any extra data (cloned objects, trivial computations, etc.) in Ruby.
* Do _NOT_ resolve VALUE into constituent elements unless they reach the function where the elements are needed or it is absolutely necessary. Passing around a VALUE in the C/C++ core is much more convienient than passing around `void*` pointers which point to an array of matrix elements.

Basically, follow a practice of 'once you enter C, never look back!'.

If you have something more in mind, discuss it in the issue tracker or
on
[this](https://groups.google.com/forum/#!topic/sciruby-dev/OJxhrGG309o)
thread.

## C/C++ style guide

This section is a work in progress.

* Use camel_case notation for arguments. No upper case.
* Write a brief description of the arguments that your function
  receives in the comments directly above the function.
* Explicitly state in the comments any anomalies that your function
  might have. For example, that it does not work with a certain
  storage or data type.

## Documentation

There are two ways in which NMatrix is being documented: guides and
comments, which are converted with RDoc into the documentation seen in
[sciruby.com](http://sciruby.com).

If you want to write a guide on how to use NMatrix to solve some
problem or simply showing how to use one of its features, write it as
a wiki page and send an e-mail on the [mailing list][1]. We're working
to improve this process.

If you aren't familiar with RDoc syntax,
[this is the official documentation](http://docs.seattlerb.org/rdoc/RDoc/Markup.html).

## Making new nmatrix extensions

From version 0.2, NMatrix supports extensions, all of which can be
hosted from the main nmatrix repo.

Refer to
[this blog post ](http://wlevine.github.io/2015/06/15/releasing-multiple-gems-with-c-extensions-from-the-same-repository.html)
to see how to do that in case you want to write your own extension for
nmatrix.

## Conclusion

This guide was heavily based on the
[Contributing to Ruby on Rails guide](http://edgeguides.rubyonrails.org/contributing_to_ruby_on_rails.html).

[1]: https://groups.google.com/forum/?fromgroups#!forum/sciruby-dev
[2]: https://github.com/sciruby/nmatrix/issues?sort=created&state=open
