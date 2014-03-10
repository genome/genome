-- Revert instrument_data_attribute

BEGIN;

DROP TABLE IF EXISTS instrument.data_attribute;

COMMIT;
