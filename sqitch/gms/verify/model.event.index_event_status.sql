-- Verify model.event.index_event_status

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'event_status_index';

ROLLBACK;
