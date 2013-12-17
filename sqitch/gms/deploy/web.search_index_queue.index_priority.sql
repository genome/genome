-- Deploy web.search_index_queue.priority
-- requires: web_search_index_queue

BEGIN;

CREATE INDEX siq_priority on web.search_index_queue using btree (priority);

COMMIT;
