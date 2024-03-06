# (C) Copyright 1996-2019 ECMWF.
# (C) Copyright 2024 KNMI.

configure_file(env_vars               env_vars                 COPYONLY)
configure_file(namelists              namelists                COPYONLY)

# === Specific for m7 tests =============================================

configure_file(m7/m7_aero.nml         namelist_m7_aero         COPYONLY)
configure_file(m7/m7_aerowvldiag.nml  namelist_m7_aerowvldiag  COPYONLY)
configure_file(m7/m7_chem.nml         namelist_m7_chem         COPYONLY)
configure_file(m7/m7_compoemis.nml    namelist_m7_compoemis    COPYONLY)

# === Parts that are replicated from test_compo_fc  ====================

configure_file(aer/aer_yaerodesc.nml  namelist_aer_yaerodesc   COPYONLY)
configure_file(aer/aer_chem.nml       namelist_aer_chem        COPYONLY)
configure_file(aer/aer_compoemis.nml  namelist_aer_compoemis   COPYONLY)
