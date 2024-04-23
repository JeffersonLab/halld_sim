# Default makefile for gxtwist
# this just invokes make in the gelhad directory and
# then make with Makefile.bms in this directory.

DIRS += gelhad hdds

.PHONY: all install clean pristine relink env current_dir

.PHONY: $(DIRS)

all install clean pristine relink env: $(DIRS) current_dir

$(DIRS):
	if test ! -d gelhad ; then ln -s ../HDGeant/gelhad . ; fi
	$(MAKE) -C $@ $(MAKECMDGOALS)

current_dir:
	$(MAKE) -- -lHDDM
	$(MAKE) -- -lHDGEOMETRY
	$(MAKE) -f Makefile.bms $(MAKECMDGOALS)

-lHDDM:
	touch -- $@
	true

-lHDGEOMETRY:
	touch -- $@
	true

clean:
	rm -rf paw.metafile flukaerr.dat geant.hbook gxtwist.rz gelhad hdgeant.hddm
	rm -f hddsGeant3_tagger.F hdds/hddsGeant4_tagger.F
	rm -f -- -lHDDM -lHDGEOMETRY

#
# DO NOT DELETE
