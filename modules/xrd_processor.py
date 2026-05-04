"""modules/xrd_processor.py — XRD module entry point for app.py."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from xrd import MODULE_INFO, run, parse_xrd_file, COMMON_WAVELENGTHS
from xrd.cod_api import (search_by_elements, search_by_name, search_by_formula,
                          fetch_cif, get_stick_pattern, SORT_OPTIONS)
from xrd.mp_api  import (search_by_elements  as mp_search_elements,
                          search_by_formula   as mp_search_formula,
                          search_by_name      as mp_search_name,
                          fetch_cif           as mp_fetch_cif,
                          validate_api_key    as mp_validate_key)
from xrd.cif_cache import get_cache, cached_fetch_cod, cached_fetch_mp
