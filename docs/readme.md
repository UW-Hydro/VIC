# VIC documentation README

This directory contains that builds the VIC website (http://vic.readthedocs.io). Below is a bit of information on how to edit, test, and serve the documentation.

## Style
The VIC documentation source files are written in [Markdown](https://help.github.com/articles/markdown-basics/), and configured with a single YAML configuration file (`mkdocs.yml`). Mardown also supports some html so that is an option if pure markdown can't get the job done.

## Requirements

To edit the documentation, all you need is a text editor.

To build the documentation, there are two requirements:
- [Python](https://www.python.org/) version 2.6, 2.7, 3.3 or 3.4.
- [`mkdocs`](http://www.mkdocs.org/) project to build its documentation.

## Building and serving the documentation locally

After editing the VIC documentation, and before committing it to the git repository, you'll want to build and serve the docs on your local machine.

#### build
from the top level of the VIC repository, run:

`mkdocs build`

#### serve
from the top level of the VIC repository, run:

`mkdocs serve`

For more information on how to interact with the built docs, check out the `mkdocs` [documentation](http://www.mkdocs.org/#getting-started).

## Read The Docs

The VIC documentation is served by [Read the Docs](https://readthedocs.io/) at http://vic.readthedocs.io. This allows us to provide multiple versions of the VIC documentation, served simultaneously from the same location.  

Currently, there are three builds (versions) scheduled for the VIC documentation:

1.  [`master`](http://vic.readthedocs.io/en/master/) - this represents the docs on VIC's `master` branch
1.  [`develop`](http://vic.readthedocs.io/en/develop/) - this represents the docs on VIC's `develop` branch
1.  `VIC.${tag}` - this represents the docs for individual tags in VIC's history.
