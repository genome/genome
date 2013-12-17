-- Verify subject_role

BEGIN;

SELECT id, name
FROM subject.role
WHERE FALSE;

ROLLBACK;
