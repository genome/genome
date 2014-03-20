-- Deploy instrument.fragment_library_unique_name
-- requires: instrument_fragment_library

BEGIN;

ALTER TABLE instrument.fragment_library ADD CONSTRAINT fragment_library_unique_full_name UNIQUE (full_name);

COMMIT;
