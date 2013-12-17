-- Verify subject_project_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'subject.project', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'subject.project', 'SELECT')::int;

ROLLBACK;
