

vmeclib :
	@cd ../VMEC2000; rm -f Release/*.o Release/*.mod Release/*.a; $(MAKE) lib_release

iwrap : 
	iwrap -f VMEC_IMAS_ACTOR.yaml

xml : 
	@xmllint --noout indata.xsd indata.xml
	@xmllint --noout --schema indata.xsd indata.xml

help :
	@echo "Simplifies building various codes for IMAS"
