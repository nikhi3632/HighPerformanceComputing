target 			:= main
objects			:= main.o functions.o

CC				:= icc
CFLAGS			:= -Wall -Werror -diag-disable=10441

RM 				:= rm -rf

all				: $(target)

depedencies 	:= $(patsubst %.o, %.d, $(objects)) 
-include $(depedencies)
DEPENDENCYFLAGS := -MMD -MF $(@:.o=.d)

main: $(objects)
	$(CC) $(CFLAGS) -o $@ $^

%o: %.c
	$(CC) $(CFLAGS) -c $< $(DEPENDENCYFLAGS)

clean:
	$(RM) $(target) $(objects) $(depedencies)
