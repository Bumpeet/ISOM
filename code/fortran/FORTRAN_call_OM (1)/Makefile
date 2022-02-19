# Dependency files
# Add them here or they will be cleaned
DEPENDENCIES = \
External_res.txt \
Makefile \
*.csv \
*.mos \
*.mo \
*.f \

CLEAN = `ls | grep -w -v -f deps.tmp`

.PHONY : test clean getdeps failingtest

test:
	@echo
	@echo Running tests...
	@$(TEST) $(TESTFILES)

# Cleans all files that are not listed as dependencies
clean:
	@echo $(DEPENDENCIES) | sed 's/ /\\\|/g' > deps.tmp
	@rm -rf $(CLEAN)

getdeps:
	@echo $(DEPENDENCIES) | sed 's/ /\\\|/g' > deps.tmp
	@echo $(CLEAN) | sed -r 's/deps.txt|deps.tmp//g' | sed 's/ / \\\n/g' > deps.txt
	@echo Dependency list saved in deps.txt.
	@echo Copy the list from deps.txt and add it to the Makefile @DEPENDENCIES

failingtest:
	@echo
	@echo Running failing tests...
	@echo
	@$(TEST) $(FAILINGTESTFILES)