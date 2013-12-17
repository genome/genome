-- Revert web.search_index_queue.index_priority

BEGIN;

DROP INDEX web.siq_priority;

COMMIT;
