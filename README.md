# obiwan Documentation

Please visit [obiwan on locally built Docs](https://legacysurvey.github.io/obiwan)

The Docs are built locally from the instructions below. The galsim install is so complicated that after getting .travis.yml working I didn't want to repeat for the Sphinx/Read the Docs build script.

This is taken from Luca Sbardella and Ryan Dale's instructions:
http://lucasbardella.com/blog/2010/02/hosting-your-sphinx-docs-in-github
https://daler.github.io/sphinxdoc-test/includeme.html

# Instructions

## Setup

1. Given a repo "github.com/legacysurvey/obiwan" containing all the Sphinx files to do "make html", the following will add a branch "gh-pages" to the repo that only contains the Sphinx build products you create on your local computer with "make html", etc.
2. git clone the repo to "./obiwan", change this line in the Sphinx Makefile to "BUILDDIR = ../../obiwandoc", commit and push to remote
3. make a directory "./obiwandocs", cd there, git clone the repo naming it as "html", cd to the html directory. e.g. you are now here "./obiwandocs/html"
4. do these
 * git branch gh-pages
 * git symbolic-ref HEAD refs/heads/gh-pages
 * rm .git/index
 * git clean -fdx
 * touch .nojekyll
5. "cd ../../obiwan/docs && make html". All the Sphinx build outputs should be written to "./../obiwandocs/"
6. "cd ../../obiwandocs/html && git add . && git commit -m 'first commit to gh-pages' && git push origin gh-pages"
7. Your docs are now live here: legacysurvey.github.com/obiwan

## How to update the docs

1. "git clone https://github.com/legacysurvey/obiwan.git && mkdir obiwandocs && cd obiwandocs && git clone https://github.com/legacysurvey/obiwan.git html && cd html && git checkout gh-pages"
2. cd obiwan/docs && make html
3. cd ../../obiwandocs/html && git add -u :/ && git commit -m "" && git push origin gh-pages"

License
=======

obiwan is free software licensed under a 3-clause BSD-style license. For details see
the ``LICENSE.rst`` file.
