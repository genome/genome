-- Verify subject_user_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'subject."user"', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'subject."user"', 'SELECT')::int;

ROLLBACK;
