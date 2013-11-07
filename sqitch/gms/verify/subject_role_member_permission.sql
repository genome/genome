-- Verify subject_role_member_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'subject.role_member', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'subject.role_member', 'SELECT')::int;

ROLLBACK;
