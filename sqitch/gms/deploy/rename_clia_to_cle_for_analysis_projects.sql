-- Deploy rename_clia_to_cle_for_analysis_projects
-- requires: add_clia_fields_to_analysis_projects

BEGIN;

ALTER TABLE config.analysis_project ADD COLUMN is_cle BOOLEAN NOT NULL DEFAULT FALSE;

COMMIT;
