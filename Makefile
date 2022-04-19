FLAGS = -Wall

all: jamming_it.c random.c
	gcc random.c $(FLAGS) -c
	gcc jamming_it.c random.o $(FLAGS) -o JAMMING_IT -lm -lgsl -lgslcblas
