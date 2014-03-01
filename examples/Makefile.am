EXAMPLES = simple_test simple_higher_rank_test 4photons 6photons ttbarh

THREADEXAMPLES = thread_6photons thread_4photons

EXTRA_PROGRAMS =  $(EXAMPLES) $(THREADEXAMPLES)

EXTRA_DIST = mynum.frm mynumhr.frm 4photons.frm 6photons.frm ttbarh.frm \
	simple_test.sh simple_higher_rank_test.sh 4photons.sh 6photons.sh \
	ttbarh.sh \
	simple_test.py simple_higher_rank_test.py 4photons.py 6photons.py \
	ttbarh.py

simple_test_SOURCES = simple_test.cc mynum.cc
simple_higher_rank_test_SOURCES = simple_higher_rank_test.cc mynumhr.cc
4photons_SOURCES = 4photons.cc 4photons_num.cc 4photons_init.cc
6photons_SOURCES = 6photons.cc 6photons_num.cc 6photons_init.cc
thread_6photons_SOURCES = thread_6photons.cc 6photons_num.cc 6photons_init.cc
thread_4photons_SOURCES = thread_4photons.cc 4photons_num.cc 4photons_init.cc
ttbarh_SOURCES = ttbarh.cc ttbarh_num.cc ttbarh_init.cc

noinst_HEADERS = mynum.hh mynumhr.hh 4photons_num.hh 6photons_num.hh \
	ttbarh_num.hh

if EXNINJASTATIC
LIBNINJA = $(libdir)/libninja.a
else
LIBNINJA = -L$(libdir) -lninja
endif

OTHEREXAMPLESCPPFLAGS = -I$(includedir)
OTHEREXAMPLESLDADD = $(LIBNINJA)  $(AVHFLAGS)  $(LTFLAGS)  $(QUADRUPLEFLAGS)

simple_test_LDADD = $(OTHEREXAMPLESLDADD)
simple_test_CPPFLAGS = $(OTHEREXAMPLESCPPFLAGS)

simple_higher_rank_test_LDADD = $(OTHEREXAMPLESLDADD)
simple_higher_rank_test_CPPFLAGS = $(OTHEREXAMPLESCPPFLAGS)

4photons_LDADD = $(OTHEREXAMPLESLDADD)
4photons_CPPFLAGS = $(OTHEREXAMPLESCPPFLAGS)

6photons_LDADD = $(OTHEREXAMPLESLDADD)
6photons_CPPFLAGS = $(OTHEREXAMPLESCPPFLAGS)

thread_6photons_LDADD = $(OTHEREXAMPLESLDADD) $(PTHREAD_CFLAGS) $(PTHREAD_LIBS)
thread_6photons_CXXFLAGS = $(OTHEREXAMPLESCPPFLAGS) $(PTHREAD_CFLAGS)

thread_4photons_LDADD = $(OTHEREXAMPLESLDADD) $(PTHREAD_CFLAGS) $(PTHREAD_LIBS)
thread_4photons_CXXFLAGS = $(OTHEREXAMPLESCPPFLAGS) $(PTHREAD_CFLAGS)

ttbarh_LDADD = $(OTHEREXAMPLESLDADD)
ttbarh_CPPFLAGS = $(OTHEREXAMPLESCPPFLAGS)

.PHONY: examples
examples: $(EXAMPLES)

.PHONY: thread-examples
thread-examples: $(THREADEXAMPLES)
