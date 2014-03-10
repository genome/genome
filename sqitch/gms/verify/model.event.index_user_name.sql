-- Verify model.event.index_user_name

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'event_user_name_index';

ROLLBACK;
