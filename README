NTL  -- a library for doing numbery theory --  version 11.6.0
Release date: 2025.11.07

Author: Victor Shoup (victor@shoup.net)

NTL is open-source software distributed under the terms of the GNU Lesser
General Public License (LGPL) version 2.1 or later.  See the file
doc/copying.txt for complete details on the licensing of NTL.

Documentation is available in the file doc/tour.html, which can be viewed with
a web browser.

For a detailed guide to installation, please see the appropriate documentation: 
   * doc/tour-unix.html for unix systems
   * doc/tour-win.html for Windows systems

The latest version of NTL is available at http://www.shoup.net.

---------------------------------
QUICK START GUIDE FOR UNIX
---------------------------------

In the src directory, you can use the basic pattern of 
   ./configure [ variable=value ]...
   make 
   make install
Note that NTL's configure script does not follow the normal autotools syntax 
for passing options to the script.
Each option is of the form variable=value.
The most common and important options:

   GMP_PREFIX=/path/to/gmp/installation  
   - by defaut, NTL relies on GMP, the GNU multi-precision library, and you will 
     typically have to specify where GMP is installed
   - the header file gmp.h should be located in $(GMP_PREFIX)/include, 
     and the library file(s) libgmp.{a,so} should be in $(GMP_PREFIX)/lib
   - note that for best performance, you shoul build GMP from source,
     rather than rely on a system installation, which is often not tuned
     to your particular machine -- especially on common x86 Linux distributions,
     where GMP is often compiled to target a generic x86 instruction set

   PREFIX=/path/to/ntl/installation
   - typically, you will have to specify where you want to install NTL
   - "make install" will copy NTL header files to $(PREFIX)/include/NTL, 
     library file(s) libntl.{a,so} to $(PREFIX)/lib, documentation files to
     $(PREFIX)/doc/NTL, and a pkg-config file ntl.pc to $(PREFIX)/lib/pkgconfig 
   
   DEF_PREFIX=/default/installation/path
   - as a short cut, you can set DEF_PREFIX so that GMP_PREFIX and PREFIX will 
     default to that value
   - for example, if you are installing both GMP and NTL to $HOME/.local,
     just set DEF_PREFIX=$HOME/.local
   
   SHARED=on
   - this defaults to "off", but you should set it to "on" if you want
     to build NTL as a shared (dynamic) library (in addition to a static 
     library)


See ../doc/config.txt for all the options.

Before installing NTL:
----------------------

You can run
   make check
to make sure everything works. This can take a few minutes.

After installing NTL:
---------------------

You can copy $(PREFIX)/include/NTL/USER_MAKEFILE.txt to "makefile" in your
current working directory.  Then you can compile the program "foo.cpp" into the
excutable "foo" by executing "make foo". 
- The makefile will compile foo.cpp with all the right flags and compiler 
  options to ensure a compilation that is
  consistent with NTL's build, and to ensure that dependencies (header 
  files and libraries) are located (at compile, link, and run time).  
- It's useful to look at this makefile to see what compiler flags are used 
  and to read the comments in the file to see why they are used on your
  particular system.  
- Of course, you can use this makefile as a starting point for more complex 
  projects.

You can look at $(PREFIX)/include/NTL/CONFIG_LOG.txt for a summary of all the
configuration options that were used in this build of NTL.

You can peruse the documentation in $(PREFIX)/doc.  Open
$(PREFIX)/doc/tour.html to get the main page of the documentation.

You can use the pkg-config file $(PREFIX)/lib/pkgconfig/ntl.pc to intregerate 
NTL into a build system such as Cmake/Meson. 


    
