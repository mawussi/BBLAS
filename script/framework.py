#! /usr/bin/env python
# -*- coding: utf-8 -*-

##
#
# @file framework.py
#
# @brief Framework used to control the installation of BBLAS and dependencies.
#
# BBLAS is a software package provided by Univ. of Manchester,
# Univ. of Tennessee.
#
# @version 1.0.0
# @author Julie Langou
# @author Mathieu Faverge
# @author Samuel D. Relton
# @date 2016-04-14
#
#
# University of Tennessee ICL License
#
# -- Innovative Computing Laboratory
# -- Electrical Engineering and Computer Science Department
# -- University of Tennessee
# -- (C) Copyright 2008-2010
#
# Redistribution  and  use  in  source and binary forms, with or without
# modification,  are  permitted  provided  that the following conditions
# are met:
#
# * Redistributions  of  source  code  must  retain  the above copyright
#   notice,  this  list  of  conditions  and  the  following  disclaimer.
# * Redistributions  in  binary  form must reproduce the above copyright
#   notice,  this list of conditions and the following disclaimer in the
#   documentation  and/or other materials provided with the distribution.
# * Neither  the  name of the University of Tennessee, Knoxville nor the
#   names of its contributors may be used to endorse or promote products
#   derived from this software without specific prior written permission.
#
# THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
# LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
##

from utils import runShellCommand, writefile, killfiles, fixpaths
import sys
import os
import getopt
import string

class Framework:
    """ This class takes care of the BBLAS installation. """

    # set default values
    ##
    # Installation directory
    ##
    prefix       = os.getcwd()+"/install"               # The install directory
    ##
    # Build directory
    ##
    build        = os.getcwd()+"/build"                 # The build directory
    ##
    # Make command
    ##
    make         = "make -j"                    # the "make" command
    ##
    # Platform
    ##
    plat         = "LINUX"
    ##
    # Fortran mangling
    ##
    mangling     = ""                        # sets the mangling type for Fortran to C calls
    ##
    # Whether or not to use MAGMA
    ##
    usemagma     = 0
    ##
    # Whether or not to use CUDA
    ##
    usecuda      = 0                         # whether or not to compile with CuBLAS
    ##
    # Whether or not to use MKL
    ##
    usemkl       = 0                         # whether or not to compile with MKL
    ##
    # Whether to use a proxy for downloading
    ##
    proxy        = 0
    ##
    # Command used for downloading
    ##
    downcmd      = ""                        # the command used to download stuff
    ##
    # Do we need TMG
    ##
    needtmg      = 0
    ##
    # Do we need to download BLAS
    ##
    downblas     = 0                         # whether or not to download reference blas
    ##
    # Do we need to download CBLAS
    ##
    downcblas    = 0                         # whether or not to download reference Cblas
    ##
    # Do we need to download LAPACK
    ##
    downlapack   = 0                         # whether or not to download reference lapack
    ##
    # Do we need to download LAPACKE
    ##
    downlapc     = 0                         # whether or not to download reference lapack c interface
    ##
    # Do we need to download TMG
    ##
    downtmg      = 0                         # whether or not to download reference tmglib
    ##
    # Do we need to download lintest
    ##
    downlintest  = 1                         # whether or not to download reference lintest
    ##
    # Do we run tests after compilation
    ##
    testing      = 1                         # whether on not to compile and run BBLAS test programs
    ##
    # Do we build the doxygen documentation
    ##
    documentation = 1
    ##
    # Version of lapack
    ##
    lapversion   = "lapack-3.6.0"
    ##
    # ranlib command
    ##
    ranlib       = ""                        # the ranlib command
    ##
    # Do we clean up afterwards
    ##
    clean        = 0
    ##
    # Do we use src option, to create makefile with no compilation
    ##
    src          = 0
    ##
    # Directory of installer
    ##
    installerdir = ""
    ##
    # Provide verbose information
    ##
    verbose      = 0
    ##
    # Is LAPACK installed
    ##
    lapackinstalled = 0;

    ##
    # Start up framework to parse arguments etc.
    ##
    def __init__(self, argv, config):

        print "="*40
        print "Setting up the framework\n"

        ##
        # Configuration file holding parameters
        ##
        self.config = config

        # parse input arguments
        self.parse_args(argv)

        if self.clean:
            self.cleanup()
            sys.exit()

        if self.config.cc == '':
            #  if no C compiler is provided
            #  raise an error
            print """
A C compiler is required to compile the Plasma library.
You can specify it using the --cc flag."""
            sys.exit()

        if self.config.fc == "":
            #  if no Fortran compiler is provided
            #  raise an error
            print """
A Fortran compiler is required to compile the core_lapack library.
You can specify it using the --fc flag."""
            sys.exit()

        if self.prefix == "" :
            ##
            # Where to install
            ##
            self.prefix = "./"
        self.prefix = fixpaths(self.prefix)
        if(not os.path.isdir(self.prefix)):
            print"Creating directory", self.prefix
            os.mkdir(self.prefix)
        print 'Install directory is...', self.prefix

        if self.build == "":
            ##
            # Where to build everything
            ##
            self.build = "./"
        self.build = fixpaths(self.build)
        if(not os.path.isdir(self.build)):
            print"Creating directory",self.build
            os.mkdir(self.build)
        print 'Build directory is...', self.build

        ##
        # Directory of installation
        ##
        self.installerdir = os.getcwd()
        os.chdir(self.build)

        # CUSTOM CHECKS
        self.check_cc()
        self.check_fc()
        self.set_mangling()
        self.set_download()
        self.set_ranlib()
        self.detect_compilers()
        self.check_linking()

        print """
**************************************
*        Compilation Options
**************************************
"""
        print 'The C compiler is... ', self.config.cc
        print 'C flags are... ', self.config.ccflags
        print 'The Fortran compiler is... ', self.config.fc
        print 'Fortran flags are... ', self.config.fcflags
        print 'Ar flags are... ', self.config.arflags

        if (self.downblas > 1) :
            print 'BLAS library is... REFERENCE BLAS ( To download )'
        else:
            print 'BLAS library is... ', self.config.blaslib

        if (self.downcblas == 2) :
            print 'CBLAS library is... REFERENCE CBLAS ( To download )'
        elif (self.config.cblaslib == ""):
            if (self.downcblas == 1):
                print 'CBLAS library is... Check if included in Blas Library and download if it is not'
            else:
                print 'CBLAS library is... Check if included in Blas Library'
        else :
            print 'CBLAS library is... ', self.config.cblaslib

        if (self.downlapack == 2) :
            print 'LAPACK library is... REFERENCE LAPACK ( To download )'
        elif (self.config.lapacklib == ""):
            if (self.downlapack == 1):
                print 'LAPACK library is... Check if included in Blas Library and download if it is not'
            else:
                print 'LAPACK library is... Check if included in Blas library'
        else:
            self.downlapack = -1
            print 'LAPACK library is... ', self.config.lapacklib

        if (self.downlapc == 2) :
            print 'LAPACK C Interface library is... REFERENCE LAPACKE INTERFACE ( To download )'
        elif (self.config.lapclib == ""):
            if (self.downlapc == 1):
                print 'LAPACKE Interface is... Check if included in Blas Library and download if it is not'
            else:
                print 'LAPACKE Interface library is... Check if included in Lapack library'
        else:
            print 'LAPACK C Interface library is... ', self.config.lapclib

        if self.config.cudadir == '':
            print """
*************************************************************************
*  CUDA library not specified. Will continue compilation without CUDA.  *
*************************************************************************
"""
        if self.config.magmadir == '':
            print """
*************************************************************************
* MAGMA library not specified. Will continue compilation without MAGMA. *
*************************************************************************
"""
        return

    ##
    # Display help information for how to use the installer.
    ##
    def usage(self):
          print "="*40
          print """
   BBLAS configuration script version %d.%d.%d
   The script will try to figure out some of the features of your system
   and the location of few other tools required for
   the installation of the BBLAS package and the other software
   packages that it requires.
   It is strongly recommended that you help the script by providing info
   through the flags listed below.

   -h or --help        : display this help and exit

   --prefix=[DIR]      : install files in DIR [./install]

   --build=[DIR]       : libraries are built in DIR [./build]
                         Contains log, downloads and builds.

   --cc=[CMD]          : the C compiler. [gcc]

   --fc=[CMD]          : the Fortran compiler. [gfortran]

   --cflags=[FLAGS]    : the flags for the C compiler [-02]

   --fflags=[FLAGS]    : the flags for the Fortran compiler [-O2]

   --ldflags_c=[flags] : loader flags when main program is in C. Some
                         compilers (e.g. PGI) require different
                         options when linking C main programs to
                         Fortran subroutines and vice-versa

   --ldflags_fc=[flags]: loader flags when main program is in
                         Fortran. Some compilers (e.g. PGI) require
                         different options when linking Fortran main
                         programs to C subroutines and vice-versa.
                         If not set, ldflags_fc = ldflags_c.

   --make=[CMD]        : the make command [make]

   --blaslib=[LIB]     : a BLAS library
                         (path should be absolute if --prefix is used)

   --cblaslib=[LIB]    : a CBLAS library
                         (path should be absolute if --prefix is used)

   --lapacklib=[LIB]   : a Lapack library
                         (path should be absolute if --prefix is used)

   --lapclib=[LIB]     : path to a LAPACK C interface.
                          (path should be absolute if --prefix is used)

   --cudadir=[DIR]     : path to CUDA installation folder e.g. "/usr/local/cuda".
                          (If not used then software is compiled without CUDA).

   --magmadir=[DIR]    : path to MAGMA installation folder e.g. "/usr/local/magma".
                          (If not used then software is compiled without MAGMA).

   --downblas          : Download and install reference BLAS.

   --downcblas         : Download and install reference CBLAS.

   --downlapack        : Download and install reference LAPACK.

   --downlapc          : Download and install reference LAPACK C Interface.

   --downall           : Download and install all missing external libraries.
                         If you don't have access to wget or no network
                         connection, you can provide the following packages
                         in the directory build/download:
        http://netlib.org/blas/blas.tgz
        http://www.netlib.org/blas/blast-forum/cblas.tgz
        http://www.netlib.org/lapack/%s.tgz

   --[no]testing       : enables/disables the testings. All externals
                         libraries are required and tested if enabled.
                         Enabled by default.

   --[no]documentation : enables/disables the generation of documentation.
                         Requires doxygen and graphviz to be installed.
                         Enabled by default.

   --disable-f90       : to disable the compilation of the f90 interface.

   --clean             : cleans up the installer directory.

   --src               : Generates a make.inc for BBLAS developers with
                         options given. If some external libraries are
                         not available, they are automatically downloaded
                         and installed in the prefix directory.
                         Testing step is also deactivated by default.
   """ % (self.config.version[0], self.config.version[1], self.config.version[2], self.lapversion)
          print "="*40
          sys.exit()

    ##
    # Parse command line arguments and set them appropriately.
    ##
    def parse_args(self, argv):
        """ Parse input argument to get compilers etc. from command line. """

        if len(argv) == 1:
            self.usage()
            sys.exit(0)

        try:
          opts, args = getopt.getopt(argv[1:], "hvp:b:n:",
                                         ["help", "prefix=", "build=",
                                          "cc=", "fc=", "cflags=", "fflags=",
                                          "ldflags_c=", "ldflags_fc=", "arflags=", "make=",
                                          "blaslib=", "cblaslib=", "lapacklib=", "lapclib=",
                                          "cudadir=", "magmadir=",
                                          "noopt=", "downblas", "downcblas", "downlapack",
                                          "downlapc", "downall", "verbose", "disable-f90",
                                          "testing", "notesting",
                                          "documentation", "nodocumentation", "clean", "src"])

        except getopt.error, msg:
          print msg
          print "for help use --help"
          sys.exit(2)

        if len(args) > 0 :
            print 'Too many arguments : ', args
            print "for help use --help"
            sys.exit(2);

        # process options
        for o, a in opts:
            if o in ("-h", "--help"):
                self.usage()
                sys.exit(0)
            else:
                if o == '--clean':
                    self.clean = 1
                    return
                elif o in ('-p', '--prefix'):
                    self.prefix = a
                elif o in ('-b', '--build'):
                    self.build = a
                elif o == '--cflags':
                    self.config.ccflags = a
                elif o=='--fflags':
                    self.config.fcflags = a
                elif o=='--noopt':
                    self.config.noopt = a
                    print 'NOOPT flags are ', a
                elif o=='--make':
                    self.make = a
                elif o=='--cc':
                    self.config.cc = a
                elif o=='--fc':
                    self.config.fc = a
                elif o == '--blaslib':
                    self.config.blaslib = fixpaths(a)
                    if "mkl" in o:
                        self.usemkl = 1
                elif o == '--downblas':
                    self.downblas = 2
                elif o == '--cblaslib':
                    self.config.cblaslib = fixpaths(a)
                elif o == '--downcblas':
                    self.downcblas = 2
                elif o == '--lapacklib':
                    self.config.lapacklib = fixpaths(a)
                elif o=='--cudadir':
                    self.config.cudadir = fixpaths(a)
                    self.usecuda = 1
                elif o=='--magmadir':
                    self.config.magmadir = fixpaths(a)
                    self.usemagma = 1
                elif o == '--downlapack':
                    self.downlapack = 2
                elif o == '--lapclib':
                    self.config.lapclib = fixpaths(a)
                elif o == '--downlapc':
                    self.downlapc = 2
                elif o == '--testing':
                    self.testing = 1
                elif o == '--notesting':
                    self.testing = 0
                elif o == '--documentation':
                    self.documentation = 1
                elif o == '--nodocumentation':
                    self.documentation = 0
                elif o == '--ldflags_c':
                    self.config.ldflags_c = a
                elif o == '--ldflags_fc':
                    self.config.ldflags_fc = a
                elif o == '--disable_f90':
                    self.config.withf90 = 0
                elif o == '--downall':
                    self.downblas   = max(1, self.downblas  )
                    self.downcblas  = max(1, self.downcblas )
                    self.downlapack = max(1, self.downlapack)
                    self.downlapc   = max(1, self.downlapc  )
                    self.downtmg    = max(1, self.downtmg   )
                elif o == '--src':
                    self.downblas   = max(1, self.downblas  )
                    self.downcblas  = max(1, self.downcblas )
                    self.downlapack = max(1, self.downlapack)
                    self.downlapc   = max(1, self.downlapc  )
                    self.downtmg    = max(1, self.downtmg   )
                    self.testing = 0
                    self.src = 1
                elif o == '--arflags':
                    self.config.arflags = a
                elif (o in ('-v', '--verbose')):
                    self.verbose = 1
                else :
                    print "Unknown option : ", o
                    sys.exit()
        # Set correctly downloads
        if (((self.config.blaslib == "") and (self.downblas > 0))
            or (self.config.blaslib == "download") ):
            self.config.blasdir = self.prefix
            self.downblas = max(1, self.downblas)
        else :
            self.downblas = 0

        if (((self.config.cblaslib == "") and (self.downcblas > 0))
            or (self.config.cblaslib == "download" )):
            self.config.cblasdir = self.prefix
            self.downcblas = max(1, self.downcblas)
        else :
            self.downcblas = 0

        if (((self.config.lapacklib == "") and (self.downlapack > 0))
            or (self.config.lapacklib == "download" )):
            self.config.lapackdir = self.prefix
            self.downlapack = max(1, self.downlapack)
        else :
            self.downlapack = 0

        if (((self.config.lapclib == "") and (self.downlapc > 0))
            or (self.config.lapclib == "download" )):
            self.config.lapcdir = self.prefix
            self.downlapc = max(1, self.downlapc)
        else :
            self.downlapc = 0

        if (self.config.ldflags_fc == "") and (self.config.ldflags_c):
            self.config.ldflags_fc = self.config.ldflags_c

    ##
    # Check if C compiler works
    ##
    def check_cc(self):
        """ checks if cc works """
        # simply generates a C program containing a couple of calls
        # to MPI routines and checks if the compilation and execution
        # are succesful
        print 'Checking if cc works...',
        sys.stdout.flush()
        # generate
        writefile('tmpc.c',"""
            #include <stdio.h>
            int main(int argc, char **argv){
            int iam;
            fprintf(stdout, \"success\" );fflush(stdout);
            return 0;
            }\n""")

        # compile
        ccomm = self.config.cc+" "+self.config.ccflags+" "+self.config.ldflags_c+" -o tmpc "+os.path.join(os.getcwd(),"tmpc.c")
        (output, error, retz) = runShellCommand(ccomm)

        if retz:
            print '\n\nCOMMON: C compiler not working! aborting...'
            print 'stderr:\n','*'*40,'\n',error,'\n','*'*40
            sys.exit()

        # run
        comm = './tmpc'
        (output, error, retz) = runShellCommand(comm)
        if retz:
            print '\n\nCOMMON: cc not working! aborting...'
            print 'error is:\n','*'*40,'\n',error,'\n','*'*40
            sys.exit()

        # cleanup
        killfiles(['tmpc.c','tmpc'])
        print 'yes'
        return 0;

    ##
    # Check if Fortran compiler works
    ##
    def check_fc(self):
        """ check if the Fortran compiler works """
        # simply generates a F77 program and checks if the compilation and execution
        # are succesful
        print "Checking if the Fortran compiler works...",
        sys.stdout.flush()
        # generate
        writefile("tmpf.f","""
      program ftest
      integer i
      print*,'success'
      stop
      end\n""")

        # compile
        fcomm = self.config.fc+' '+self.config.fcflags+" "+self.config.ldflags_fc+' -o tmpf '+'tmpf.f'
        (output, error, retz) = runShellCommand(fcomm)

        if retz:
            print '\n\nCOMMON: the Fortran compiler is not working! aborting...'
            print 'error is:\n','*'*40,'\n',error,'\n','*'*40
            sys.exit()

        # run
        comm = './tmpf'
        (output, error, retz) = runShellCommand(comm)
        if retz:
            print '\n\nCOMMON: the Fortran compiler is not working! aborting...'
            print 'error is:\n','*'*40,'\n',error,'\n','*'*40
            sys.exit()

        # cleanup
        killfiles(['tmpf.f','tmpf','tmpf.o'])
        print 'yes'

        return 0;

    ##
    # Set mangler for Fortran
    ##
    def set_mangling(self):
        """ Sets the INTFACE variable in Bmake.inc """
        # This one generates a program equivalent to that in BLACS/INSTALL
        # that checks the mangling in FORTRAN function symbols
        print 'Setting Fortran mangling...',
        sys.stdout.flush()
        writefile('tmpf.f',"""
      program intface
      external c_intface
      integer i
      call c_intface(i)
      stop
      end\n""")
        writefile('tmpc.c',"""
      #include <stdio.h>
      void c_intface_(int *i){fprintf(stdout, \"-DADD_\");fflush(stdout);}
      void c_intface(int *i){fprintf(stdout, \"-DNOCHANGE\");fflush(stdout);}
      void c_intface__(int *i){fprintf(stdout, \"-DfcIsF2C\");fflush(stdout);}
      void C_INTFACE(int *i){fprintf(stdout, \"-DUPCASE\");fflush(stdout);}\n""")

        ccomm = self.config.cc+' '+self.config.ccflags+' -c tmpc.c -o tmpc.o'
        fcomm = self.config.fc+' '+self.config.fcflags+' '+self.config.ldflags_fc+' tmpf.f tmpc.o -o xintface'

        (output, error, retz) = runShellCommand(ccomm)
        if retz:
            print '\n\nCOMMON: in set_mangling: cannot compile'
            print 'error is:\n','*'*40,'\n',error,'\n','*'*40
            sys.exit()

        (output, error, retz) = runShellCommand(fcomm)
        if retz:
            print '\n\nCOMMON: in set_mangling: cannot compile'
            print 'error is:\n','*'*40,'\n',error,'\n','*'*40
            sys.exit()

        comm = os.path.join(os.getcwd(),'xintface')
        (output, error, retz) = runShellCommand(comm)
        if retz:
            print '\n\nCOMMON: in set_mangling: cannot run xintface'
            print 'error is:\n','*'*40,'\n',error,'\n','*'*40
            sys.exit()

        self.mangling = output
        killfiles(['xintface', 'tmpf.f', 'tmpf.o', 'tmpc.c', 'tmpc.o'])

        print self.mangling
        return 1;


    ##
    # Check if C main can be linked to Fortran subroutine
    ##
    def check_linking(self):
        """ Check if C main can be linked to Fortran subroutine """

        # This one checks if the linking command works out of the box or
        # if any specific flag is required. For example if the linker if the
        # Intel FORTRAN compiler, then the "-nofor_main" is usually required.
        # This function only checks if linker works but does not automatically
        # detect the required flags
        print 'Checking loader...',
        sys.stdout.flush()
        writefile('tmpf.f',"""
      subroutine fsub()
      write(*,*)'success'
      stop
      end\n""")
        writefile('tmpc.c',"""
      #if defined ADD_
      #define fsub fsub_
      #elif defined NOCHANGE
      #define fsub fsub
      #elif defined fcIsF2C
      #define fsub fsub_
      #elif defined UPCASE
      #define fsub FSUB
      #endif
      void main(){
      fsub();}\n""")

        ccomm = self.config.cc+' '+self.config.ccflags+' '+self.mangling+' -c -o tmpc.o tmpc.c'
        fcomm = self.config.fc+' '+self.config.fcflags+' -c -o tmpf.o tmpf.f'
        lcomm = self.config.fc+' '+self.config.ldflags_fc+' '+self.config.ld_fcmain+' -o lnk tmpf.o tmpc.o'

        (output, error, retz) = runShellCommand(ccomm)
        if retz:
            print '\n\nCOMMON: in check_linking: cannot compile'
            print 'command is: ',ccomm
            print 'error is:\n','*'*40,'\n',error,'\n','*'*40
            sys.exit()

        (output, error, retz) = runShellCommand(fcomm)
        if retz:
            print '\n\nCOMMON: in check_linking: cannot compile'
            print 'command is: ',fcomm
            print 'error is:\n','*'*40,'\n',error,'\n','*'*40
            sys.exit()

        (output, error, retz) = runShellCommand(lcomm)
        if retz:
            print """\n\nCOMMON: in check_linking: cannot link
            Cannot link a C main program to a Fortran77 subroutine
            Make sure that the appropriate flags are passed to the linker."""
            print 'command is: ',lcomm
            print 'error is:\n','*'*40,'\n',error,'\n','*'*40
            sys.exit()


        killfiles(['lnk', 'tmpf.f', 'tmpf.o', 'tmpc.c', 'tmpc.o'])

        print 'works'
        return 1;



    ##
    # Figure out how to download files
    ##
    def set_download(self):
        """ Figures out how to download files """
        print 'Setting download command...'
        wget = 0
        urllib = 0
        # JULIE : Cut proxy stuff...was causing problems (see scalapack installer if you want it back)
        if urllib == 0:
            # if urllib2 is not present checks if wget is present
            # in the PATH and if yes it sets the download command
            # to be wget
            print "Checking availablility of wget...",
            path=str(os.getenv('PATH')).split(os.pathsep)
            for i in path:
                if (os.path.isfile(os.path.join(i,'wget'))):
                    print "available"
                    wget = 1
                    break
            if wget:
                # test wget
                print "Testing wget...",
                comm = 'wget --tries=2 --timeout=5 http://www.netlib.org/lapack/index'
                (output, error, retz) = runShellCommand(comm)
                if(retz != 0):
                    print 'not working.'
                    wget = -1
                else:
                    print "working"
                    self.downcmd="wget"
                    os.remove("index")
                    return
            else:
                # wget not available
                print "not available"
                wget=0


    ##
    # Sets the ranlib command
    ##
    def set_ranlib(self):
        """ Sets the ranlib command """
        # Some systems don't have the ranlib command (e.g. SGIs).
        # In the case where ranlib is not present in the PATH,
        # echo is used instead of ranlib
        print "Setting ranlib command...",

        path=str(os.getenv('PATH')).split(os.pathsep)
        for i in path:
            if os.path.isfile(os.path.join(i,'ranlib')):
                self.config.ranlib=os.path.join(i,'ranlib')
                print self.config.ranlib
                return

        for i in path:
            if os.path.isfile(os.path.join(i,'echo')):
                self.config.ranlib=os.path.join(i,'echo')
                print self.config.ranlib
                return

    ##
    # Detect compilers
    ##
    def detect_compilers(self):
        """ Tries to detect the compilers type """
        # By users experience it is known which compiler flags are required
        # in some cases. This function tries to detect which compilers are used
        # and sets the flags accordingly

        print 'Detecting Fortran compiler...',
        if self.fc_is_intel():
            # The Intel FORTRAN compiler requires -nofor_main flag
            # for the linking and the -mp flag to maintain the
            # floating-point precision
            self.config.fcflags    += ' -diag-disable vec -fltconsistency -fp_port'
            self.config.ldflags_c  += ' '   # used to link
            self.config.ldflags_fc += ' '
            self.config.ld_fcmain   = ' -nofor_main'
            self.config.noopt      += ' -mp'
            self.testing = 0; # Cannot compile lintest with fc_main option
            print 'Intel'
        elif self.fc_is_gnu():
            print 'GNU'
            self.config.ld_fcmain   = ''
        elif self.fc_is_xlf():
            self.config.fcflags    += ' -qstrict -qthreaded'
            self.config.ld_fcmain   = ''
            print 'IBM'
        elif self.fc_is_pgi():
            self.config.ldflags_c  += ''
            self.config.ldflags_fc += ''
            self.config.ld_fcmain   = ' -Mnomain'
            self.testing = 0; # Cannot compile lintest with fc_main option
        else:
            self.config.compiler = "Unknown"
            print 'unknown'

        print 'Detecting C compiler...',
        if self.cc_is_intel():
            self.config.compiler = "Intel"
            self.config.ccflags    += ' -diag-disable vec'
            print 'Intel'
        elif self.cc_is_gnu():
            self.config.compiler = "GNU"
            print 'GNU'
        elif self.cc_is_xlc():
            self.config.compiler = "XLC"
            self.config.ccflags    += ' -qstrict -qthreaded'
            print 'IBM'
        elif self.cc_is_pgi():
            self.config.compiler = "PGI"
            print 'PGI'
        else:
            print 'unknown'

        print 'Selected C compiler flags: '+self.config.ccflags
        print 'Selected Fortran compiler flags: '+self.config.fcflags
        print 'Selected loader flags (C main): '+self.config.ldflags_c
        print 'Selected loader flags (Fortran main): '+self.config.ldflags_fc
        return


    ##
    # Check if Fortran compiler is from Intel
    ##
    def fc_is_intel(self):
        comm = self.config.fc+' -V'
        (output, error, retz) = runShellCommand(comm)
        isifort = string.find(error,'Intel(R) Fortran')
        if isifort != -1:
            return 1

        return 0

    ##
    # Check if Fortran compiler is GNU
    ##
    def fc_is_gnu(self):
        comm = self.config.fc+' -v'
        (output, error, retz) = runShellCommand(comm)
        isifort = string.find(error,'gcc')
        if isifort != -1:
            return 1

        return 0

    ##
    # Check if Fortran compiler is XLF
    ##
    def fc_is_xlf(self):
        comm = self.config.fc
        (output, error, retz) = runShellCommand(comm)
        isifort = string.find(output,'xlf')
        if isifort != -1:
            return 1

        return 0

    ##
    # Check if Fortran compiler is PGI
    ##
    def fc_is_pgi(self):
        comm = self.config.fc+' -V'
        (output, error, retz) = runShellCommand(comm)
        isifort = string.find(error,'pgifc')
        if isifort != -1:
            return 1
        isifort = string.find(error,'pgif95')
        if isifort != -1:
            return 1
        isifort = string.find(error,'Portland')
        if isifort != -1:
            return 1
        isifort = string.find(output,'pgifc')
        if isifort != -1:
            return 1
        isifort = string.find(output,'pgif95')
        if isifort != -1:
            return 1
        isifort = string.find(output,'Portland')
        if isifort != -1:
            return 1

        return 0

    ##
    # Check C compiler is Intel
    ##
    def cc_is_intel(self):
        comm = self.config.cc+' -V'
        (output, error, retz) = runShellCommand(comm)
        isifort = string.find(error,'Intel(R) C')
        if isifort != -1:
            return 1

        return 0

    ##
    # Check C compiler is GNU
    ##
    def cc_is_gnu(self):
        comm = self.config.cc+' -v'
        (output, error, retz) = runShellCommand(comm)
        isifort = string.find(error,'gcc')
        if isifort != -1:
            return 1

        return 0

    ##
    # Check C compiler is XLC
    ##
    def cc_is_xlc(self):
        comm = self.config.cc+' -qversion'
        (output, error, retz) = runShellCommand(comm)
        isifort = string.find(output,'XL')
        if isifort != -1:
            return 1

        return 0

    ##
    # Check C compiler is PGI
    ##
    def cc_is_pgi(self):
        comm = self.config.cc+' -V'
        (output, error, retz) = runShellCommand(comm)
        isifort = string.find(error,'pgicc')
        if isifort != -1:
            return 1
        isifort = string.find(error,'Portland')
        if isifort != -1:
            return 1
        isifort = string.find(output,'pgicc')
        if isifort != -1:
            return 1
        isifort = string.find(output,'Portland')
        if isifort != -1:
            return 1

        return 0



    ##
    # Continue with installation
    ##
    def resume(self):

        if self.src :
            print """
******************************************************
******************************************************

Installation of all librairies required by BBLAS is complete.

You can copy the following file in your src directory to compile
BBLAS :\n"""+os.path.join(self.build,'make.inc')+"""\n

WARNING: the prefix directory has been set to same directory than
external libraries.  If you want to install BBLAS in another
directory, don't forget to change it.

We strongly recommend the use of hwloc library for a better binding of
the threads.  It will be automatically detected if pkg-config is set
to find it.

Log messages are in the\n"""+os.path.join(self.build,'log')+"""
directory.

The\n"""+self.build+"""/downloads
directory contains all tarballs of downloaded librairies.

The\n"""+self.build+"""
directory contains the source code of the libraries that have been
installed. It can be removed once you have copy the make.inc file.

******************************************************
******************************************************
"""
        else:
            print """
******************************************************
******************************************************

BBLAS installation completed.

Your BLAS library is               : """+self.config.blasdir+"""\n
Your CBLAS library is              : """+self.config.cblasdir+"""\n
Your LAPACK library is             : """+self.config.lapackdir+"""\n
Your LAPACK C Interface library is : """+self.config.lapcdir+"""\n

You can add this line to your environment configuration file to
use pkg-config with plasma :
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:"""+self.prefix+"""/lib/pkgconfig

Log messages are in the\n"""+os.path.join(self.build,'log')+"""
directory.

The BBLAS testing programs are in:\n"""+os.path.join(self.build,'plasma/testing')+"""

The\n"""+self.build+"""
directory contains the source code of the libraries
that have been installed. It can be removed at this time.
"""
        return 0

    ##
    # Clean up installer directory
    ##
    def cleanup(self):
        """ Cleans up the installer directory """

        print "Cleaning up...",
        sys.stdout.flush()

        builddir = os.path.join(self.build)

        comm = 'rm -rf '+builddir
        #+' '+libdir+' '+logdir
        (output, error, retz) = runShellCommand(comm)

        print "done."
