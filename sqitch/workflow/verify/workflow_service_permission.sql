-- Verify workflow_service_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'workflow.service', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'workflow.service', 'SELECT')::int;

ROLLBACK;
