###############################################################################
# Sloan Digital Sky Survey III (SDSS-III)
# Code for product: cosmoxi2d
#
# This Makefile & all Makefiles in this product are GNU make compliant.
# Please help keep them that way.  See
# http://www.gnu.org/software/make/manual/make.html
#
# $Id: Makefile 138170 2012-09-11 19:25:55Z weaver $
#
###############################################################################
GSL_DIR ?= /usr/local
MYINCLUDE = -I$(GSL_DIR)/include
#
# Depending on the system, this might have to be lib64
#
MYLIB = -L$(GSL_DIR)/lib
#MYLIB = -L$(GSL_DIR)/lib64
LDFLAGS = $(MYLIB) -lm -lgsl -lgslcblas
CC = gcc
LD = gcc
CFLAGS = $(MYINCLUDE)
BIN = $(COSMOXI2D_DIR)/bin
OBJS = misc.o besselroots.o peakbackground.o matsubaraQRandSPT.o \
	linearpk.o ptraw.o pktoxi.o pktoxilin.o pktoxiLPT.o pktovin.o \
	pktosig.o xistream.o xistreamLPT.o xicombine.o cosmoxi2d.o
#
# SDSS-III Makefiles should always define this target.
#
all: $(BIN)/cosmoxi2d
#
# Individual objects
#
linearpk.o: misc.o
ptraw.o: linearpk.o matsubaraQRandSPT.o
pktoxilin.o: linearpk.o besselroots.o
pktoxiLPT.o: misc.o ptraw.o besselroots.o
pktoxi.o: linearpk.o ptraw.o besselroots.o
pktovin.o: linearpk.o ptraw.o besselroots.o
pktosig.o: linearpk.o ptraw.o besselroots.o pktovin.o
xistream.o: misc.o pktoxi.o pktovin.o pktosig.o
xistreamLPT.o: misc.o pktoxiLPT.o
xicombine.o: xistreamLPT.o xistream.o
cosmoxi2d.o: misc.o peakbackground.o linearpk.o ptraw.o pktoxi.o \
	pktoxilin.o pktoxiLPT.o pktovin.o pktosig.o xistream.o xistreamLPT.o xicombine.o
#
# Executable
#
$(BIN)/cosmoxi2d: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)
#
# GNU make pre-defines $(RM).  The - in front of $(RM) causes make to
# ignore any errors produced by $(RM).
#
clean:
	- $(RM) *.o *~ $(BIN)/cosmoxi2d
