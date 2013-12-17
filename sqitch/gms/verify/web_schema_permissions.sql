-- Verify web_schema_permissions

BEGIN;

select 1/has_schema_privilege('genome', 'web', 'CREATE')::int;
select 1/has_schema_privilege('gms-user', 'web', 'USAGE')::int;

ROLLBACK;
