-- Revert web_search_index_queue

BEGIN;

DROP TABLE IF EXISTS web.search_index_queue;

COMMIT;
