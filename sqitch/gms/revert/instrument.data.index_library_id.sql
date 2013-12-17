-- Revert instrument.data.index_library_id

BEGIN;

DROP INDEX instrument.instrument_data_library_id_index;

COMMIT;
