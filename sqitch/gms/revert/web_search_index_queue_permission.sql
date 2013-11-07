-- Revert web_search_index_queue_permission

BEGIN;

REVOKE ALL ON TABLE web.search_index_queue FROM "gms-user";

COMMIT;
