-- Verify workflow_instance_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'workflow.instance', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'workflow.instance', 'SELECT')::int;

ROLLBACK;
