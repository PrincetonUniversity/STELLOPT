

vmeclib :
	@cd ../VMEC2000; rm -f Release/*.o Release/*.mod Release/*.a; $(MAKE) lib_release

actor_vmec : 
	iwrap -f VMEC_IMAS_ACTOR.yaml

xml_vmec : 
	@xmllint --noout indata.xsd indata.xml
	@xmllint --noout --schema indata.xsd indata.xml

help :
	@echo "Simplifies building various codes for IMAS"
	@echo "  vmeclib: Clean build of vmec and library"
	@echo "  xml_vmec: Checks VMEC XML files"
	@echo "  actor_vmec: Builds VMEC actor using iwrap"
