CXX = 		g++ -g -fmessage-length=0 -w

SOURCE =	src

OUTPUT =	bin

LD_DIRS = 	-lpthread

DCP_LIBS =	-IOGDF -LOGDF/_release -lOGDF

all:	decomp0 decomp1
	
decomp0:
	$(CXX) -o $(OUTPUT)/decomp0 $(SOURCE)/decomp0.cpp $(DCP_LIBS) $(LD_DIRS)
	
decomp1:
	$(CXX) -o $(OUTPUT)/decomp1 $(SOURCE)/decomp1.cpp $(LD_DIRS) $(DCP_LIBS)

clean:
	rm -f $(SOURCE)/*.o $(OUTPUT)/*
