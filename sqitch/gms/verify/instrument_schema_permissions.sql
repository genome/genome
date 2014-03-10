-- Verify instrument_schema_permissions

BEGIN;

select 1/has_schema_privilege('genome', 'instrument', 'CREATE')::int;
select 1/has_schema_privilege('gms-user', 'instrument', 'USAGE')::int;

ROLLBACK;
