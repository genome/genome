-- Deploy web_search_index_queue_permission
-- requires: web_search_index_queue

BEGIN;

REVOKE ALL ON TABLE web.search_index_queue FROM PUBLIC;
REVOKE ALL ON TABLE web.search_index_queue FROM genome;
GRANT ALL ON TABLE web.search_index_queue TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE web.search_index_queue TO "gms-user";

COMMIT;
