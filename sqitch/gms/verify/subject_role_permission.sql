-- Verify subject_role_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'subject.role', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'subject.role', 'SELECT')::int;

ROLLBACK;
