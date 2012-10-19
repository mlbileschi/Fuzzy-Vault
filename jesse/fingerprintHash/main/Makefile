CXXFLAGS = -O5 #may throw illegal instruction
LDFLAGS = -I ../common/NTL-mac/include -L ../common/NTL-mac/lib
LINKFLAGS = -lntl
SRCDIR = ../common/src
INCDIR = ../common/inc


all:  ProbDistroGA.o IOUtils.o evaluatePerformance.o WelchBerlekamp.o VaultMethod.o vault.o VaultTriangles.o VaultBF.o runTest.o unitTests.o unitTestsVault.o unitTestsGA.o VaultLetsTrySomething.o main

#for windows
#clean: 
#	`cmd` //C del *.o

#for unix/mac
clean:
	rm *.o *~ RUNME ../common/inc/*~ ../common/src/*~

ProbDistroGA.o: $(SRCDIR)/ProbDistroGA.cpp $(INCDIR)/ProbDistroGA.h $(INCDIR)/runTest.h
	g++ $(CXXFLAGS) -c $(SRCDIR)/ProbDistroGA.cpp $(LDFLAGS)

WelchBerlekamp.o: $(SRCDIR)/WelchBerlekamp.cpp $(INCDIR)/WelchBerlekamp.h
	g++ $(CXXFLAGS) -c $(SRCDIR)/WelchBerlekamp.cpp $(LDFLAGS)

VaultMethod.o: $(SRCDIR)/VaultMethod.cpp $(INCDIR)/VaultMethod.h $(INCDIR)/common.h
	g++ $(CXXFLAGS) -c $(SRCDIR)/VaultMethod.cpp $(LDFLAGS)

VaultTriangles.o:  $(SRCDIR)/VaultTriangles.cpp $(INCDIR)/VaultMethod.h
	g++ $(CXXFLAGS) -c $(SRCDIR)/VaultTriangles.cpp $(LDFLAGS)

VaultBF.o:  $(SRCDIR)/VaultBF.cpp $(INCDIR)/VaultMethod.h
	g++ $(CXXFLAGS) -c $(SRCDIR)/VaultBF.cpp $(LDFLAGS)

vault.o: $(SRCDIR)/vault.cpp $(INCDIR)/vault.h $(INCDIR)/common.h $(SRCDIR)/VaultMethod.cpp $(INCDIR)/VaultMethod.h $(SRCDIR)/VaultTriangles.cpp
	g++ $(CXXFLAGS) -c $(SRCDIR)/vault.cpp $(LDFLAGS)

IOUtils.o:  $(SRCDIR)/IOUtils.cpp $(INCDIR)/IOUtils.h $(INCDIR)/common.h $(INCDIR)/runTest.h
	g++ $(CXXFLAGS) -c $(SRCDIR)/IOUtils.cpp $(LDFLAGS)

evaluatePerformance.o: $(SRCDIR)/evaluatePerformance.cpp $(INCDIR)/evaluatePerformance.h $(SRCDIR)/IOUtils.cpp
	g++ $(CXXFLAGS) -c $(SRCDIR)/evaluatePerformance.cpp $(LDFLAGS)

runTest.o: $(INCDIR)/runTest.h $(SRCDIR)/runTest.cpp $(INCDIR)/evaluatePerformance.h $(INCDIR)/IOUtils.h 
	g++ $(CXXFLAGS) -c $(SRCDIR)/runTest.cpp $(LDFLAGS)

unitTests.o: $(INCDIR)/common.h $(SRCDIR)/unitTests.cpp $(INCDIR)/unitTests.h
	g++ $(CXXFLAGS) -c $(SRCDIR)/unitTests.cpp $(LDFLAGS)

unitTestsVault.o: $(INCDIR)/common.h $(SRCDIR)/unitTestsVault.cpp $(INCDIR)/unitTests.h
	g++ $(CXXFLAGS) -c $(SRCDIR)/unitTestsVault.cpp $(LDFLAGS)

unitTestsGA.o: $(INCDIR)/common.h $(SRCDIR)/unitTestsGA.cpp $(INCDIR)/unitTests.h
	g++ $(CXXFLAGS) -c $(SRCDIR)/unitTestsGA.cpp $(LDFLAGS)

VaultLetsTrySomething.o: $(SRCDIR)/VaultLetsTrySomething.cpp $(INCDIR)/VaultMethod.h
	g++ $(CXXFLAGS) -c $(SRCDIR)/VaultLetsTrySomething.cpp $(LDFLAGS)

main.o: main.cpp $(INCDIR)/evaluatePerformance.h $(INCDIR)/IOUtils.h $(INCDIR)/unitTests.h
	g++ $(CXXFLAGS) -c main.cpp $(LDFLAGS)

main: main.o 
	g++ $(CXXFLAGS) IOUtils.o ProbDistroGA.o main.o evaluatePerformance.o vault.o VaultMethod.o WelchBerlekamp.o VaultTriangles.o VaultBF.o runTest.o unitTests.o unitTestsVault.o unitTestsGA.o VaultLetsTrySomething.o -o RUNME $(LDFLAGS) $(LINKFLAGS)

