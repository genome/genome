-- Revert subject.project_part.index_project_id_label

BEGIN;

DROP INDEX subject.project_part_project_label_index;

COMMIT;
