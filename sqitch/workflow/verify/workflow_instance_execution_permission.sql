-- Verify workflow_instance_execution_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'workflow.instance_execution', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'workflow.instance_execution', 'SELECT')::int;

ROLLBACK;
