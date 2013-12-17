-- Verify web_task_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'web.task', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'web.task', 'SELECT')::int;

ROLLBACK;
