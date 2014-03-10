-- Verify instrument.data_attribute.index_instrument_data_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'instrument_data_id_index';

ROLLBACK;
