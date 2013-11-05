-- Deploy subject.project_part.label
-- requires: subject_project_part

BEGIN;

CREATE INDEX idx_s_pp_l on subject.project_part using btree (label);

COMMIT;
