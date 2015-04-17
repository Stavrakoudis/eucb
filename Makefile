PROGRAMNAME=eucb
INSTDIR  = ~/bin
CXX      = c++
CFLAGS   = -m64 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE
CXXFLAGS = -m64
LDFLAGS  = -m64

all:
	(cd src; make alldist)

clean:
	(cd src; make cleandist)
	(cd bin; rm -f $(PROGRAMNAME))

install: all
	(cp -f bin/$(PROGRAMNAME) $(INSTDIR) )
