-- Revert subject.project.index_name

BEGIN;

DROP INDEX subject.project_name_index;

COMMIT;
