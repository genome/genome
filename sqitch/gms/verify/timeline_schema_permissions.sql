-- Verify timeline_schema_permissions

BEGIN;

select 1/has_schema_privilege('genome', 'timeline', 'CREATE')::int;
select 1/has_schema_privilege('gms-user', 'timeline', 'USAGE')::int;

ROLLBACK;
