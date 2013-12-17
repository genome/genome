-- Verify disk_schema_permissions

BEGIN;

select 1/has_schema_privilege('genome', 'disk', 'CREATE')::int;
select 1/has_schema_privilege('gms-user', 'disk', 'USAGE')::int;

ROLLBACK;
