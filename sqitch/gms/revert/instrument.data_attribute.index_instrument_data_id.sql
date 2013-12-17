-- Revert instrument.data_attribute.index_instrument_data_id

BEGIN;

DROP INDEX instrument.instrument_data_id_index;

COMMIT;
