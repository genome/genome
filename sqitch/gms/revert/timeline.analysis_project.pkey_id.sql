-- Revert timeline.analysis_project.pkey_id.sql

BEGIN;
    ALTER TABLE timeline.analysis_project DROP CONSTRAINT analysis_project_pkey;
COMMIT;
