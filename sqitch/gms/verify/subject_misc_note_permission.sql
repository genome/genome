-- Verify subject_misc_note_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'subject.misc_note', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'subject.misc_note', 'SELECT')::int;

ROLLBACK;
