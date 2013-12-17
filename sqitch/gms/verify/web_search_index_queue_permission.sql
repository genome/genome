-- Verify web_search_index_queue_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'web.search_index_queue', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'web.search_index_queue', 'SELECT')::int;

ROLLBACK;
