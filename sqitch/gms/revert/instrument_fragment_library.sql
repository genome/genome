-- Revert instrument_fragment_library

BEGIN;

DROP TABLE IF EXISTS instrument.fragment_library;

COMMIT;
