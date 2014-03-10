-- Deploy config.analysis_project.index_status
-- requires: analysis_project

BEGIN;

CREATE INDEX analysis_project_status_idx ON config.analysis_project USING btree (status);

COMMIT;
