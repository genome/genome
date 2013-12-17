-- Revert subject.project_part.index_project_id_role

BEGIN;

DROP INDEX subject.project_part_project_role_index;

COMMIT;
