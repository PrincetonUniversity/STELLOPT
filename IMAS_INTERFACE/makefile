.PHONEY : vmeclib actor_vmec xml_vmec clean_vmec beams3dlib actor_beams3d xml_beams3d clean_beams3d help

#####  VMEC

vmeclib :
	@cd ../VMEC2000; $(MAKE) lib_release

actor_vmec : vmeclib
	iwrap -f VMEC_IMAS_ACTOR.yaml -i $(ACTOR_FOLDER)

xml_vmec : 
	@xmllint --noout indata.xsd indata.xml
	@xmllint --noout --schema indata.xsd indata.xml

clean_vmec :
	@cd ../VMEC2000; rm -f Release/*.o Release/*.mod Release/*.a

#####  BEAMS3D

beams3dlib :
	@cd ../BEAMS3D; $(MAKE) lib_release

actor_beams3d : beams3dlib
	iwrap -f BEAMS3D_IMAS_ACTOR.yaml -i $(ACTOR_FOLDER)

xml_beams3d : 
	@xmllint --noout beams3d_input.xsd beams3d_input.xml
	@xmllint --noout --schema beams3d_input.xsd beams3d_input.xml

clean_beams3d :
	@cd ../BEAMS3D; rm -f Release/*.o Release/*.mod Release/*.a


#####  HELP message

help :
	@echo "Simplifies building various codes for IMAS"
	@echo "  vmeclib: Clean build of vmec and library"
	@echo "  xml_vmec: Checks VMEC XML files"
	@echo "  actor_vmec: Builds VMEC actor using iwrap"
	@echo "  beams3dlib: Clean build of vmec and library"
	@echo "  xml_beams3d: Checks VMEC XML files"
	@echo "  actor_beams3d: Builds VMEC actor using iwrap"
