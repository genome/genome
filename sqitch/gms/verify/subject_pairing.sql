-- Verify subject_pairing

BEGIN;

SELECT id, created_at, updated_at, created_by, control_subject_id,
    experimental_subject_id, analysis_project_id
FROM subject.pairing
WHERE FALSE;

ROLLBACK;
