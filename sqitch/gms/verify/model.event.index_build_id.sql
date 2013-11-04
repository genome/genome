-- Verify model.event.index_build_id

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'event_build_id_index';

ROLLBACK;
