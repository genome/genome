-- Deploy subject.project.name
-- requires: subject_project

BEGIN;

CREATE INDEX project_name_index on subject.project using btree (name);

COMMIT;
