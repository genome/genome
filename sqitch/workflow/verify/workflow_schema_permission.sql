-- Verify workflow_schema_permission

BEGIN;

select 1/has_schema_privilege('genome', 'workflow', 'CREATE')::int;
select 1/has_schema_privilege('gms-user', 'workflow', 'USAGE')::int;

ROLLBACK;
