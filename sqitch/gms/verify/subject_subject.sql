-- Verify subject_subject

BEGIN;

SELECT subject_id, subclass_name, name
FROM subject.subject
WHERE FALSE;

ROLLBACK;
