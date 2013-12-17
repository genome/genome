-- Revert subject.project_part.index_part_id

BEGIN;

DROP INDEX subject.idx_subject_project_part_part_id;

COMMIT;
