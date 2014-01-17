-- Revert add_clia_fields_to_analysis_projects

BEGIN;

  ALTER TABLE config.analysis_project DROP COLUMN is_cle;
  ALTER TABLE config.analysis_project DROP COLUMN run_as;

COMMIT;
