.PHONEY : vmeclib actor_vmec xml_vmec clean_vmec help

vmeclib :
	@cd ../VMEC2000; $(MAKE) lib_release

actor_vmec : vmeclib
	iwrap -f VMEC_IMAS_ACTOR.yaml -i $(ACTOR_FOLDER)

xml_vmec : 
	@xmllint --noout indata.xsd indata.xml
	@xmllint --noout --schema indata.xsd indata.xml

clean_vmec :
	@cd ../VMEC2000; rm -f Release/*.o Release/*.mod Release/*.a

help :
	@echo "Simplifies building various codes for IMAS"
	@echo "  vmeclib: Clean build of vmec and library"
	@echo "  xml_vmec: Checks VMEC XML files"
	@echo "  actor_vmec: Builds VMEC actor using iwrap"
