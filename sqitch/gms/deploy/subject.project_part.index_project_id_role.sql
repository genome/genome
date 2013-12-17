-- Deploy subject.project_part.project_id_role
-- requires: subject_project_part

BEGIN;

CREATE INDEX project_part_project_role_index on subject.project_part using btree (project_id, role);

COMMIT;
