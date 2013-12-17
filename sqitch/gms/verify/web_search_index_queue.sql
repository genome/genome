-- Verify web_search_index_queue

BEGIN;

SELECT id, subject_id, subject_class, "timestamp", priority
FROM web.search_index_queue
WHERE FALSE;

ROLLBACK;
