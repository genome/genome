-- Verify process_schema_permissions

BEGIN;

select 1/has_schema_privilege('genome', 'process', 'CREATE')::int;
select 1/has_schema_privilege('gms-user', 'process', 'USAGE')::int;

ROLLBACK;
