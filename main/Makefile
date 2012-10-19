CXX = g++
CXXFLAGS = -O5 #optimization may throw illegal instruction
LDFLAGS = -I ../common/NTL-mac/include -L ../common/NTL-mac/lib
LINKFLAGS = -lntl
SRCDIR = ../common/src
INCDIR = ../common/inc
TJEASRCDIR = ../matching_tjea/src
TJEAINCDIR = ../matching_tjea/inc


all: matching.o ProbDistroGA.o IOUtils.o EvaluatePerformance.o WelchBerlekamp.o VaultMethod.o Vault.o VaultTriangles.o VaultBF.o RunTest.o UnitTests.o UnitTestsVault.o UnitTestsGA.o VaultLetsTrySomething.o fingerprint_model.o global_parameters.o Minutiae2SFeature.o polar_util.o SFV_Util.o ConvexHull.o ExtendMatch.o SFV_VerifyMatch.o NumOverlappedMinutiae.o logfile.o Main

#for windows
#clean: 
#	`cmd` //C del *.o

#for unix/mac
clean:
	rm -f *.o *~ RUNME ../common/inc/*~ ../common/src/*~ ../doc/*~

ProbDistroGA.o: $(SRCDIR)/ProbDistroGA.cpp $(INCDIR)/ProbDistroGA.h $(INCDIR)/RunTest.h
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/ProbDistroGA.cpp $(LDFLAGS)

WelchBerlekamp.o: $(SRCDIR)/WelchBerlekamp.cpp $(INCDIR)/WelchBerlekamp.h
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/WelchBerlekamp.cpp $(LDFLAGS)

VaultMethod.o: $(SRCDIR)/VaultMethod.cpp $(INCDIR)/VaultMethod.h $(INCDIR)/Common.h
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/VaultMethod.cpp $(LDFLAGS)

VaultTriangles.o:  $(SRCDIR)/VaultTriangles.cpp $(INCDIR)/VaultMethod.h
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/VaultTriangles.cpp $(LDFLAGS)

VaultBF.o:  $(SRCDIR)/VaultBF.cpp $(INCDIR)/VaultMethod.h
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/VaultBF.cpp $(LDFLAGS)

Vault.o: $(SRCDIR)/Vault.cpp $(INCDIR)/Vault.h $(INCDIR)/Common.h $(SRCDIR)/VaultMethod.cpp $(INCDIR)/VaultMethod.h $(SRCDIR)/VaultTriangles.cpp
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/Vault.cpp $(LDFLAGS)

IOUtils.o:  $(SRCDIR)/IOUtils.cpp $(INCDIR)/IOUtils.h $(INCDIR)/Common.h $(INCDIR)/RunTest.h
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/IOUtils.cpp $(LDFLAGS)

EvaluatePerformance.o: $(SRCDIR)/EvaluatePerformance.cpp $(INCDIR)/EvaluatePerformance.h $(SRCDIR)/IOUtils.cpp
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/EvaluatePerformance.cpp $(LDFLAGS)

RunTest.o: $(INCDIR)/RunTest.h $(SRCDIR)/RunTest.cpp $(INCDIR)/EvaluatePerformance.h $(INCDIR)/IOUtils.h 
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/RunTest.cpp $(LDFLAGS)

UnitTests.o: $(INCDIR)/Common.h $(SRCDIR)/UnitTests.cpp $(INCDIR)/UnitTests.h
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/UnitTests.cpp $(LDFLAGS)

UnitTestsVault.o: $(INCDIR)/Common.h $(SRCDIR)/UnitTestsVault.cpp $(INCDIR)/UnitTests.h
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/UnitTestsVault.cpp $(LDFLAGS)

UnitTestsGA.o: $(INCDIR)/Common.h $(SRCDIR)/UnitTestsGA.cpp $(INCDIR)/UnitTests.h
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/UnitTestsGA.cpp $(LDFLAGS)

VaultLetsTrySomething.o: $(SRCDIR)/VaultLetsTrySomething.cpp $(INCDIR)/VaultMethod.h
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/VaultLetsTrySomething.cpp $(LDFLAGS)


fingerprint_model.o: $(SRCDIR)/fingerprint_model.cpp $(INCDIR)/fingerprint_model.h $(INCDIR)/random_gaussian.h
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/fingerprint_model.cpp $(LDFLAGS)

matching.o: $(TJEASRCDIR)/matching.cpp $(TJEAINCDIR)/common.h $(TJEAINCDIR)/matching.h  $(TJEAINCDIR)/SFeature.h
	$(CXX) $(CXXFLAGS) -c $(TJEASRCDIR)/matching.cpp $(LDFLAGS)

global_parameters.o: $(TJEASRCDIR)/global_parameters.cpp 
	$(CXX) $(CXXFLAGS) -c $(TJEASRCDIR)/global_parameters.cpp $(LDFLAGS)

Minutiae2SFeature.o: $(TJEASRCDIR)/Minutiae2SFeature.cpp 
	$(CXX) $(CXXFLAGS) -c $(TJEASRCDIR)/Minutiae2SFeature.cpp $(LDFLAGS)

polar_util.o: $(TJEASRCDIR)/polar_util.cpp 
	$(CXX) $(CXXFLAGS) -c $(TJEASRCDIR)/polar_util.cpp $(LDFLAGS)


SFV_Util.o: $(TJEASRCDIR)/SFV_Util.cpp 
	$(CXX) $(CXXFLAGS) -c $(TJEASRCDIR)/SFV_Util.cpp $(LDFLAGS)

ConvexHull.o: $(TJEASRCDIR)/ConvexHull.cpp 
	$(CXX) $(CXXFLAGS) -c $(TJEASRCDIR)/ConvexHull.cpp $(LDFLAGS)


ExtendMatch.o: $(TJEASRCDIR)/ExtendMatch.cpp 
	$(CXX) $(CXXFLAGS) -c $(TJEASRCDIR)/ExtendMatch.cpp $(LDFLAGS)


SFV_VerifyMatch.o: $(TJEASRCDIR)/SFV_VerifyMatch.cpp 
	$(CXX) $(CXXFLAGS) -c $(TJEASRCDIR)/SFV_VerifyMatch.cpp $(LDFLAGS)

NumOverlappedMinutiae.o: $(TJEASRCDIR)/NumOverlappedMinutiae.cpp 
	$(CXX) $(CXXFLAGS) -c $(TJEASRCDIR)/NumOverlappedMinutiae.cpp $(LDFLAGS)

logfile.o: $(TJEASRCDIR)/logfile.cpp 
	$(CXX) $(CXXFLAGS) -c $(TJEASRCDIR)/logfile.cpp $(LDFLAGS)




Main.o: Main.cpp $(INCDIR)/EvaluatePerformance.h $(INCDIR)/IOUtils.h $(INCDIR)/UnitTests.h
	$(CXX) $(CXXFLAGS) -c Main.cpp $(LDFLAGS)

Main: Main.o 
	$(CXX) $(CXXFLAGS) IOUtils.o ProbDistroGA.o Main.o EvaluatePerformance.o Vault.o VaultMethod.o WelchBerlekamp.o VaultTriangles.o VaultBF.o RunTest.o UnitTests.o UnitTestsVault.o UnitTestsGA.o VaultLetsTrySomething.o fingerprint_model.o matching.o global_parameters.o Minutiae2SFeature.o polar_util.o SFV_Util.o ConvexHull.o ExtendMatch.o SFV_VerifyMatch.o NumOverlappedMinutiae.o logfile.o -o RUNME $(LDFLAGS) $(LINKFLAGS)
	echo -e "\nexecutable file's name is RUNME\n"

