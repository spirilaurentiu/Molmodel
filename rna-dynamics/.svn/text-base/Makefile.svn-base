%.exe : %.o Makefile
	/usr/bin/c++      -Wno-deprecated   -g -headerpad_max_install_names -fPIC  $< -o $@ -L/usr/local/lib  -L$$DYLD_LIBRARY_PATH  -lSimTKsimbody  -lpthread -lm -lSimTKmolmodel -lSimTKsimbody -lSimTKmath -lSimTKcpodes -lSimTKcommon -lSimTKlapack -lpthread



#                               list of dependencies
%.o  :  %.cpp
	export DYLD_LIBRARY_PATH=$$SIMTK_PREFIX/lib
	/usr/bin/c++   -Wno-deprecated   -g -I$$SIMTK_PREFIX/include/molmodel/internal -I$$SIMTK_PREFIX/include -c  $<
	#/usr/bin/c++   -Wno-deprecated -no-cpp-precomp  -g -I/Users/samuelflores/svn/molmodel/include/molmodel/internal -I$$SIMTK_PREFIX/include  -c  $<
