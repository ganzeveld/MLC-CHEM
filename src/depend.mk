# DO NOT DELETE THIS LINE - used by make depend
messy_emdep.o: messy_main_tools.o
messy_emdep_emis.o: messy_emdep_emis_mem.o messy_main_constants_mem.o
messy_emdep_emis.o: messy_main_tools.o
messy_emdep_emis_mem.o: messy_emdep.o messy_main_constants_mem.o
messy_emdep_xtsurf.o: messy_emdep.o messy_emdep_mem.o
messy_emdep_xtsurf.o: messy_main_constants_mem.o messy_main_tools.o
messy_emdep_xtsurf_box.o: messy_emdep.o messy_emdep_emis.o
messy_emdep_xtsurf_box.o: messy_emdep_emis_mem.o messy_emdep_mem.o
messy_emdep_xtsurf_box.o: messy_emdep_xtsurf.o messy_main_constants_mem.o
messy_main_tools.o: messy_main_constants_mem.o
