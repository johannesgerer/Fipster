# SHELL := /bin/bash
LDLIBS= -lstdc++ -lm
CXX = gccfilter -a g++ $(IDIR)
CXX = g++ $(IDIR)
CC = g++
IDIR = -Isrc -I. -Iinclude
CXXFLAGS = -fdiagnostics-show-location=once \
					-fmessage-length=80 \
					-std=c++11 \
					-pipe \
					-Winvalid-pch\
					-Wall \
					#-H \
					#-Werror \



prog=src/main

precompiled=src/precompiled.h
precompiled_auto=src/precompiled_auto.h
precompile_flag=AUTOMATIC_PRECOMPILATION

tests = $(wildcard tests/*.cpp)
sources = $(wildcard src/*.cpp)

deps = $(addsuffix .d, $(sources) $(precompiled) $(tests) )

run: debug
	./$(prog)

runo: opt
	./$(prog)

debug: CXXFLAGS += -g -DDEBUG  -D_GLIBCXX_DEBUG
debug: LDFLAGS += -g
debug: LDLIBS += -ltbb_debug
debug: $(prog)
opt: CXXFLAGS += -g -O3 -DNDEBUG
opt: LDFLAGS += -g
opt: LDLIBS += -ltbb
opt: $(prog)
$(prog): LDLIBS += -lboost_system -lboost_filesystem -lboost_program_options
$(prog): $(sources:.cpp=.o)

test: unittest
	./unittest
unittest: LDLIBS += -lboost_unit_test_framework
unittest: $(tests:.cpp=.o)

cleantest:
	rm unittest -f


#forces existence of precompiled.h.gch, which in turn includes precompiled_auto.h (via included dependencies), which is generated via the precompiled_auto.h.tmp rule
$(sources:.cpp=.o): $(addsuffix .gch, $(precompiled))

#extract everything within the ifndef block and remove duplicates
#after first occurence
$(precompiled_auto).tmp:
	echo "#define $(precompile_flag)" > $@
	for i in src/*.h src/*.cpp ; do \
		sed  -e '1,/#ifndef $(precompile_flag)/ d' -e '/#endif/,$$ d' $$i; \
	done | awk ' !x[$$0]++' >>  $@
	if ! diff -N -q $@ $(precompiled_auto); then \
		mv $@ $(precompiled_auto); else rm $@; fi

$(precompiled_auto): $(precompiled_auto).tmp

%.gch: %
	$(COMPILE.cpp) -Wno-error $(OUTPUT_OPTION) $<

cleandeps:
	rm -f $(deps)

cleanobjs:
	rm -f $(prog) $(sources:.cpp=.o)

clean: cleandeps cleanobjs cleanprecomp

cleanprecomp:
	rm -f $(precompiled_auto) $(precompiled_auto).* $(precompiled).d


-include $(deps)

%.d: %
	$(COMPILE.cpp) -MG -MM '$<' \
	 -MT '$(patsubst %.h, %.h.gch, $(<:.cpp=.o)) $@' -MF '$@'
