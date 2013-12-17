-- Verify subject_schema_permissions

BEGIN;

select 1/has_schema_privilege('genome', 'subject', 'CREATE')::int;
select 1/has_schema_privilege('gms-user', 'subject', 'USAGE')::int;

ROLLBACK;
