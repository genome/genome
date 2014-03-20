-- Revert instrument.fragment_library_unique_name

BEGIN;

ALTER TABLE instrument.fragment_library DROP CONSTRAINT fragment_library_unique_full_name;

COMMIT;
