-- Verify workflow_historian_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'workflow.historian', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'workflow.historian', 'SELECT')::int;

ROLLBACK;
