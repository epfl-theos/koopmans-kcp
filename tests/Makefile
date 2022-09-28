
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

include ./ENVIRONMENT

default :
	@echo 'To run the koopmans-kcp test-suite, type at the shell prompt:'
	@echo ' '
	@echo '  make target'
	@echo ' '
	@echo 'where <target> identifies an action'
	@echo ' run-tests                   : run-tests-serial (default)'
	@echo ' run-tests-serial            : run all tests in SERIAL'
	@echo ' run-tests-parallel          : run all tests in PARALLEL (4 MPI)'
	@echo ' compare                     : compare last output with reference'
	@echo ' clean                       : clean stdout/sderr of all tests'
	@echo ' '
	@echo 'For additional advanced commands and settings please manually inspect'
	@echo 'ENVIRONMENT and Makefile files'

prolog :
	@sed "s|XXXXXX|$(ESPRESSO_ROOT)|g" < userconfig.tmp > userconfig

run-tests : run-tests-serial

run-tests-serial : run-tests-kcp-serial

run-tests-parallel : run-tests-kcp-parallel

run-tests-kcp-serial : prolog
	@if ! test -d $(ESPRESSO_TMPDIR); then mkdir $(ESPRESSO_TMPDIR); fi
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose

run-tests-kcp-parallel : prolog
	@if ! test -d $(ESPRESSO_TMPDIR); then mkdir $(ESPRESSO_TMPDIR); fi
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose

compare : compare-kcp

compare-kcp : prolog
	${TESTCODE_DIR}/bin/testcode.py --verbose compare 2>&1 | tee out.KCP.`date +%Y%m%d_%H%M%S`

benchmarks : prolog
	${TESTCODE_DIR}/bin/testcode.py --verbose make-benchmarks

clean:
# output files
	@rm -rf kcp_*/test*
# error files
	@rm -rf kcp_*/CRASH
# Cleaning KCW final files
	@rm -rf	kcp_*/out kcp_*/tmp
