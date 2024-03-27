CXX := g++
CXXFLAGS := -Wall
INCLUDES := -I.

SRCS := template_class_vector.cpp
OBJS := $(SRCS:.cpp=.o)
EXEC := template_class_vector

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