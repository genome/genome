-- Deploy remove_is_clia_from_analysis_project
-- requires: rename_clia_to_cle_for_analysis_projects

BEGIN;

ALTER TABLE config.analysis_project DROP COLUMN is_clia;

COMMIT;
