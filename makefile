VERSION:="$(shell git describe 2>/dev/null || echo '??')"

.PHONY: doc src

default: release

release:
	$(MAKE) -C src

openmpgnu:
	$(MAKE) -C src openmp=gnu

openmpintel:
	$(MAKE) -C src openmp=intel

debug:
	$(MAKE) -C src config=debug


clean:
	make -C src clean
	make -C doc clean

doc:
	make -C doc all
