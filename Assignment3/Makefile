target 			:= main 
objects			:= main.o tendency.o flux.o
shallow			:= shallow.*

CC				:= g++
CFLAGS			:= -Wall -Werror

RM 				:= rm -rf

all				: $(target)

depedencies 	:= $(patsubst %.o, %.d, $(objects)) 
-include $(depedencies)
DEPENDENCYFLAGS := -MMD -MF $(@:.o=.d)

main: $(objects)
	$(CC) $(CFLAGS) -o $@ $^

%o: %.cc
	$(CC) $(CFLAGS) -c $< $(DEPENDENCYFLAGS)

clean:
	$(RM) $(target) $(objects) $(depedencies)
	$(RM) $(shallow)
	 


