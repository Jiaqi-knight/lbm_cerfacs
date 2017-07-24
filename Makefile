# Declaration of variables
CC = g++
CC_FLAGS = -O3 -w -ansi -pedantic -fopenmp -g
LD_FLAGS = -fopenmp

COTEXT  = "\033[1;31m Compiling\033[0m\033[1m $(<F)\033[0m\n"
LITEXT  = "\033[1;31m Building \033[0m\033[1m $@\033[0m\n"

# File names
EXEC    = $(BINDIR)LEOPARD
SRCDIR = ./src/
OBJDIR = ./.Obj/
BINDIR = ./bin/
MKDIRS = $(OBJDIR) $(BINDIR)
CPP_FILES = $(wildcard $(SRCDIR)*.cpp)
OBJ_FILES = $(addprefix $(OBJDIR), $(notdir $(CPP_FILES:.cpp=.o)))

# Main target
$(EXEC): mes compil_info $(MKDIRS) $(OBJ_FILES) 
	@echo
	@printf $(LITEXT) 
	@$(CC) $(LD_FLAGS) $(OBJ_FILES) -o $(EXEC)

# To obtain object files
$(OBJDIR)%.o: $(SRCDIR)%.cpp 
	@printf $(COTEXT) 
	@$(CC) -c $(CC_FLAGS) $< -o  $@

# To remove generated files
clean: 
	@printf "\033[1;31m Cleaning object, executable and archive files\033[0m\n"
	@rm -f $(EXEC) $(OBJ_FILES) LEOPARD.tgz

# To create a tar archive of the project
tar: clean
	@printf "\033[1;31m Creating tar archive of the code \033[0m\n"
	@mkdir -p LEOPARD_student
	@cp -r Makefile doc src testcases LEOPARD_student
	@tar czf LEOPARD_student.tgz --exclude .git --exclude .Obj --exclude output LEOPARD_student
	@rm -rf LEOPARD_student

$(MKDIRS): 
	@mkdir -p $@

compil_info:
	@printf "\033[1;31m Compiler used \033[0m \033[1m    [$(CC)]\033[0m\n"
	@printf "\033[1;31m Compiling options\033[0m \033[1m [$(CC_FLAGS)]\033[0m\n" 
	@printf "\033[1;31m Linking options \033[0m \033[1m  [$(LD_FLAGS)]\033[0m\n"
	@echo

mes:
	@echo
	@echo    "       *************************************************************"
	@echo    "       *                                                           *"
	@echo    "       *     __      ______  ____    ____    ___     ____    ____  *"
	@echo    "       *    / /     / ____/ / __ \  / __ \  /   |   / __ \  / __ \ *"
	@echo    "       *   / /     / __/   / / / / / /_/ / / /| |  / /_/ / / / / / *"
	@echo    "       *  / /___  / /___  / /_/ / / ____/ / ___ | / _, _/ / /_/ /  *"
	@echo    "       * /_____/ /_____/  \____/ /_/     /_/  |_|/_/ |_| /_____/   *"
	@echo    "       *                                                           *"
	@printf "       *         \033[0m \033[1mL\033[0mattic\033[0m\033[1mE\033[0m b\033[0m\033[1mO\033[0mltzmann \033[0m \033[1mP\033[0ml\033[0m\033[1mA\033[0mtefo\033[0m\033[1mR\033[0mm \033[0m\033[1mD\033[0mevelopment         *\n"
	@echo    "       *************************************************************"
	@echo    "       *                   Copyright 2015, CERFACS                 *"
	@printf "       *************************************************************\033[0m\n"
	@echo


# To print the help message
help: mes
	@echo
	@printf "\033[1;31m Make options of LEOPARD code\033[0m\n"
	@echo
	@printf "\033[1;31m Compiler used \033[0m \033[1m    [$(CC)]\033[0m\n"
	@printf "\033[1;31m Compiling options\033[0m \033[1m [$(CC_FLAGS)]\033[0m\n" 
	@printf "\033[1;31m Linking options \033[0m \033[1m  [$(LD_FLAGS)]\033[0m\n"
	@echo	
	@printf "\033[1;31m Provided Rules: \033[0m\n"
	@printf "\033[1;31m  help         =>\033[0m\033[1m printing this help message\033[0m\n"
	@printf "\033[1;31m  clean        =>\033[0m\033[1m cleaning object, executable and archive files\033[0m\n"
	@printf "\033[1;31m  tar          =>\033[0m\033[1m creating a tar archive of the project\033[0m\n"
	@echo
