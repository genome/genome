-- Revert subject_project_part

BEGIN;

DROP TABLE IF EXISTS subject.project_part;

COMMIT;
