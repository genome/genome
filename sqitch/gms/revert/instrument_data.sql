-- Revert instrument_data

BEGIN;

DROP TABLE IF EXISTS instrument.data;

COMMIT;
