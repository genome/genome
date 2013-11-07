-- Verify web_task_params_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'web.task_params', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'web.task_params', 'SELECT')::int;

ROLLBACK;
