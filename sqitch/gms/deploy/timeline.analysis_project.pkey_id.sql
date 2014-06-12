-- Deploy timeline.analysis_project.pkey_id.sql
-- requires: analysis_project_event_logging

CREATE UNIQUE INDEX CONCURRENTLY analysis_project_pkey_idx ON timeline.analysis_project (id);
BEGIN;
    ALTER TABLE timeline.analysis_project ADD CONSTRAINT analysis_project_pkey PRIMARY KEY USING INDEX analysis_project_pkey_idx;
COMMIT;
