-- Verify subject_subject_attribute_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'subject.subject_attribute', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'subject.subject_attribute', 'SELECT')::int;

ROLLBACK;
