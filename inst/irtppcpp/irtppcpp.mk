include $(R_HOME)/etc$(R_ARCH)/Makeconf

SRCS = src/type/ghquads.cpp \
	src/utils/Input.cpp \
	src/type/dataset.cpp \
	src/estimation/estep.cpp \
	src/estimation/mstep.cpp \
	src/utils/asa111.cpp

LIBRARY = irtppcpp
SRC_DIR = src
INCLUDES = -I./$(SRC_DIR) -I./include/SPGO/include/
OBJS = $(SRCS:.cpp=.o)

all : lib$(LIBRARY).a

$(OBJS): %.o : %.h

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX1X) $(CXX1XSTD) $(CXX1XFLAGS) $(CXX1XPICFLAGS) -c $(INCLUDES) -o $@ $<

lib$(LIBRARY).a : $(OBJS)
	$(AR) rsv $@ $^

clean:
	$(RM) $(OBJS)
	$(RM) lib$(LIBRARY).a
