-- Revert subject_project

BEGIN;

DROP TABLE IF EXISTS subject.project;

COMMIT;
