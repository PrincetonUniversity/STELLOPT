# Default flags (optional)
CFLAGS += -O2 -march=native -g -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare

# Flags required by the Makefile
override CPPFLAGS += -MMD -MP

# Default target
.PHONY: all
all: clmdif.o

.PHONY: clean
clean:
	$(RM) clmdif.o clmdif.d
