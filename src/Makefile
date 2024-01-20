CC=gcc -lstdc++
CFLAGS= -Wall -Werror -Wextra -std=c++17
LIBS=-lgtest
LEAKS=leaks --atExit -- ./test

OS=$(shell uname -s)
ifeq ($(OS), Linux)
	LIBS+= -lgmock -lpthread
	LEAKS=valgrind ./test
	CC=g++
endif

all: test s21_matrix_oop.a

clean:
	rm -rf s21_matrix_oop.a *.gcda *.gcno *.info *.o
	rm -rf test_report .clang-format *.o test *.out report

test: test_build 
	./test

s21_matrix_oop.a: s21_matrix_oop.o
	ar -rcs s21_matrix_oop.a s21_matrix_oop.o

gcov_report: clean
	@$(CC) --coverage s21_matrix_oop.cpp test_matrix.cpp -o test_report $(LIBS)
	./test_report
	lcov -t "test" -o test.info --no-external -c -d .
	genhtml -o report test.info
	open report/index.html

leaks: test_build
	$(LEAKS)

clang:
	cp ../materials/linters/.clang-format .
	clang-format -n *.cpp *.h

test_build: s21_matrix_oop.o test.o
	@$(CC) $(CFLAGS) -o test s21_matrix_oop.o test_matrix.o $(LIBS)

test.o:	test_matrix.cpp
	@$(CC) $(CFLAGS) -c test_matrix.cpp $(LIBS)

s21_matrix_oop.o: s21_matrix_oop.cpp
	@$(CC) $(CFLAGS) -c s21_matrix_oop.cpp

to_clang:
	clang-format -i *.cpp *.h
