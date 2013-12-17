-- Verify workflow_plan_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'workflow.plan', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'workflow.plan', 'SELECT')::int;

ROLLBACK;
