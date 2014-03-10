-- Deploy config.analysis_project.index_name
-- requires: analysis_project

BEGIN;

CREATE INDEX analysis_project_name_idx ON config.analysis_project USING btree (name);

COMMIT;
