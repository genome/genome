-- Verify web.search_index_queue.index_priority

BEGIN;

SELECT 1/count(*) FROM pg_class WHERE relkind = 'i' and relname = 'siq_priority';

ROLLBACK;
