Mazurka Plugins
================

The Mazurka Plugins are a set of [Vamp plugins](http://www.vamp-plugins.org),
particularly for use with [Sonic Visualiser](http://www.sonic-visualiser.org).


Downloading
===========

To download with git:

```bash
git clone https://github.com/craigsapp/MazurkaPlugins
```

Or download with the zip button on the [Mazurka Plugin Github page](https://github.com/craigsapp/MazurkaPlugins).


Configuration
=============

Source code for individual plugins are found in the
[src-plugin](https://github.com/craigsapp/MazurkaPlugins/tree/master/src-plugin)
directory.  Edit the file
[Makefile.plugin](https://github.com/craigsapp/MazurkaPlugins/tree/master/Makefile.plugin#L10-L12)
to select which plugins should be included in the compiled plugin file.


Compiling
=========

The plugin Makefile is currently optimized to compile for OS X.  
To compile, type:

```bash
make
```

The first time this is run, this will download and compile two 
external libraries:

* [FFTW](http://www.fftw.org): Fast Fourier Transform library
* [Vamp SDK](http://www.vamp-plugins.org): The Vamp plugin SDK

(Note that wget needs to be installed to download.  Will update so
that curl will be instead used on OS X; otherwise, if you have 
[Homebrew](http://brew.sh) installed, type `brew install wget`
in OS X to install wget).


Installing
==========

The compiled plugin is found in `lib-osx/mazurka-plugins.dylib`.
Copy this file to `~/Library/Audio/Plugin-Ins/Vamp` and (re)start
Sonic Visualiser to used the new plugin file from the Transform menu.
To copy the plugin to the plugin directory, type:

```bash
make install
```



