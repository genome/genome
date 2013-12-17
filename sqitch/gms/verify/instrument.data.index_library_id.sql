-- Verify instrument.data.index_library_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'instrument_data_library_id_index';

ROLLBACK;
