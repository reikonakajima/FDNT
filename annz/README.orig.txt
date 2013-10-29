
ANNz Installation
=============================

The ANNz user guide is provided in postscript format in the doc/ directory.

The ANNz source code is to be found in the src/ directory. The command
'make all' issued in this directory will compile the three
executables: annz_net, annz_train and annz_test. Use 'make clean' to
remove temporary object files created during compilation.

The example directory contains a sample data set (from the Sloan
Digital Sky Survey Early Data Release), to allow you to quickly test
the functionality of ANNz. Use this with reference to the user guide.

Please feel free to contact me (Adrian Collister) at aac@ast.cam.ac.uk
with any problems etc.

Thank you for using ANNz!


--

ANNz directory structure.
--------------------------

> ls -lR annz
annz:
total 16
drwxr-xr-x  2 adrian users 4096 Nov  3 10:32 doc/
drwxr-xr-x  2 adrian users 4096 May 15 14:41 example/
-rw-r--r--  1 adrian users 2251 Nov  3 10:33 readme
drwxr-xr-x  2 adrian users 4096 Nov  3 10:35 src/

annz/doc:
total 980
-rw-r--r--  1 adrian users 996284 May 14 20:03 user_guide.ps

annz/example:
total 1492
-r--r--r--  1 adrian users 756000 May 14 20:29 sdss.ugriz.test
-r--r--r--  1 adrian users 630000 May 14 20:29 sdss.ugriz.train
-r--r--r--  1 adrian users 126000 May 14 20:29 sdss.ugriz.valid

annz/src:
total 100
-rw-r--r--  1 adrian users  2259 Nov 10 12:26 annz_net.cpp
-rw-r--r--  1 adrian users  1147 Nov 10 12:26 annz_test.cpp
-rw-r--r--  1 adrian users  1019 Nov 10 12:26 annz_train.cpp
-rw-r--r--  1 adrian users   801 Nov 10 12:26 Makefile
-rw-r--r--  1 adrian users 12193 Nov 10 12:26 network.cpp
-rw-r--r--  1 adrian users  4220 Nov 10 12:26 network.h
-rw-r--r--  1 adrian users  6813 Nov 10 12:26 test.cpp
-rw-r--r--  1 adrian users  1342 Nov 10 12:26 test.h
-rw-r--r--  1 adrian users  4115 Nov 10 12:26 test_network.cpp
-rw-r--r--  1 adrian users  1222 Nov 10 12:26 test_network.h
-rw-r--r--  1 adrian users 16804 Nov 10 12:26 train.cpp
-rw-r--r--  1 adrian users  2338 Nov 10 12:26 train.h
-rw-r--r--  1 adrian users  1287 Nov 10 12:26 train_network.cpp
-rw-r--r--  1 adrian users  1200 Nov 10 12:26 train_network.h
-rw-r--r--  1 adrian users  3347 Nov 10 12:26 util.cpp
-rw-r--r--  1 adrian users  1294 Nov 10 12:26 util.h
