-- Deploy subject.project_part.part_class_name_part_id_role
-- requires: subject_project_part

BEGIN;

CREATE INDEX project_part_part_role_index on subject.project_part using btree (part_class_name, part_id, role);

COMMIT;
