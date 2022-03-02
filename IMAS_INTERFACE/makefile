

vmeclib :
	@cd ../VMEC2000; rm -f Debug/*.o Debug/*.mod Debug/*.a; $(MAKE) lib_debug

iwrap : 
	iwrap -f VMEC_IMAS_ACTOR.yaml

xml : 
	@xmllint --noout indata.xsd indata.xml
	@xmllint --noout --schema indata.xsd indata.xml

help :
	@echo "Simplifies building various codes for IMAS"