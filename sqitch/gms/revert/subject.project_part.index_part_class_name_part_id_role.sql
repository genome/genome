-- Revert subject.project_part.index_part_class_name_part_id_role

BEGIN;

DROP INDEX subject.project_part_part_role_index;

COMMIT;
