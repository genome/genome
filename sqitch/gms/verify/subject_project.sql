-- Verify subject_project

BEGIN;

SELECT id, name
FROM subject.project
WHERE FALSE;

ROLLBACK;
