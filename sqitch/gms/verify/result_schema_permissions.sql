-- Verify result_schema_permissions

BEGIN;

select 1/has_schema_privilege('genome', 'result', 'CREATE')::int;
select 1/has_schema_privilege('gms-user', 'result', 'USAGE')::int;

ROLLBACK;
