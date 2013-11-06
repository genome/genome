-- Verify config_schema_permissions

BEGIN;

select 1/has_schema_privilege('genome', 'config', 'CREATE')::int;
select 1/has_schema_privilege('gms-user', 'config', 'USAGE')::int;

ROLLBACK;
