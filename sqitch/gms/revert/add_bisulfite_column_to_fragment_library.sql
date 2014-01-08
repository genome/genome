-- Revert add_bisulfite_column_to_fragment_library

BEGIN;

  ALTER TABLE instrument.fragement_library DROP COLUMN bisulfite_conversion;

COMMIT;
