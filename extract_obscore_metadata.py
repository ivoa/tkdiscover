"""
This script takes obscore column metadata from an installed DaCHS.
It produces source code for inclusion in tkdiscover and is only
interesting for code maintainers.
"""

from gavo import api
from gavo.utils import typeconversions

import pprint

local_columns = {"preview", "source_table"}

metadata = []
for col in api.getRD("//obscore").getById("ObsCore"):
    if col.name in local_columns:
        continue

    if col.xtype=="adql:REGION":
        # special case handled in serialisers in current: DaCHS: the
        # ugly STC-S in s_region.
        vot_type, vot_arraysize, vot_xtype = "char", "*", "adql:REGION"
    else:
        # the sane rest
        vot_type, vot_arraysize, vot_xtype = \
            typeconversions.sqltypeToVOTable(col.type)


    metadata.append({
       "name": col.name,
       "description": col.description,
       "ucd": col.ucd,
       "utype": col.utype,
       "datatype": vot_type,
       "arraysize": vot_arraysize,
       "xtype": vot_xtype,})

metadata.append(
    {'arraysize': '*',
      'datatype': 'char',
      'description': "IVOID of the originating service",
      'name': 'origin_service',
      'ucd': 'meta.ref.ivoid',
      'utype': None,
      'xtype': None})

print("OBSCORE_METADATA = \\")
pprint.pprint(metadata)
