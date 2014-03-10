-- Verify model.event.index_instrument_data_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'event_inst_data_index';

ROLLBACK;
