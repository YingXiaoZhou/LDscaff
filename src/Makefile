OBJECTS  = graph.o tools.o cal.o 
CXXFLAGS = -O3 -g -lemon -std=c++0x

matching : $(OBJECTS)
	g++ -o matching $(OBJECTS) $(CXXFLAGS)

graph.o : graph.cc
	g++ -c graph.cc $(CXXFLAGS)

tools.o : cal.cc
	g++ -c cal.cc $(CXXFLAGS)

tools.o : tools.cc
	g++ -c tools.cc $(CXXFLAGS)

.PHONY : clean

clean : 
	rm matching *.o
