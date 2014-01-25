-- Deploy add_bisulfite_column_to_fragment_library
-- requires: instrument_fragment_library

BEGIN;

  ALTER TABLE instrument.fragment_library ADD COLUMN bisulfite_conversion varchar(32);

COMMIT;
