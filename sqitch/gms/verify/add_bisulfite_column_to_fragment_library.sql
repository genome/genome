-- Verify add_bisulfite_column_to_fragment_library

BEGIN;

  SELECT bisulfite_conversion FROM instrument.fragment_library WHERE FALSE;

ROLLBACK;
