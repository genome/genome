-- Verify model_schema_permissions

BEGIN;

select 1/has_schema_privilege('genome', 'model', 'CREATE')::int;
select 1/has_schema_privilege('gms-user', 'model', 'USAGE')::int;

ROLLBACK;
