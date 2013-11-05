-- Deploy subject.project_part.part_id
-- requires: subject_project_part

BEGIN;

CREATE INDEX idx_subject_project_part_part_id on subject.project_part using btree (part_id);

COMMIT;
