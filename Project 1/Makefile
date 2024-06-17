CXX := g++
CXXFLAGS := -Wall -fopenmp -ggdb
INCLUDES := -I.

SRCS := project_1.cpp
OBJS := $(SRCS:.cpp=.o)
EXEC := project_1

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

run: 
	./$(EXEC)

clean:
	$(RM) $(EXEC) $(OBJS)

cleanimage:
	$(RM) image.png

# Second project

project:
	g++ project_2.cpp lbfgs.c -o project

runpr:
	./project

remove:
	rm project