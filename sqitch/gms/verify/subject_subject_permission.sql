-- Verify subject_subject_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'subject.subject', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'subject.subject', 'SELECT')::int;

ROLLBACK;
