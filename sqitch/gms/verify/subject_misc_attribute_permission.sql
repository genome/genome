-- Verify subject_misc_attribute_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'subject.misc_attribute', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'subject.misc_attribute', 'SELECT')::int;

ROLLBACK;
