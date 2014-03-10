-- Verify model.event.index_ref_seq_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'event_ref_seq_index';

ROLLBACK;
