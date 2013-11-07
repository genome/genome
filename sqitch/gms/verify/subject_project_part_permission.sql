-- Verify subject_project_part_permission

BEGIN;

SELECT 1/has_table_privilege('genome', 'subject.project_part', 'TRUNCATE')::int;
SELECT 1/has_table_privilege('gms-user', 'subject.project_part', 'SELECT')::int;

ROLLBACK;
