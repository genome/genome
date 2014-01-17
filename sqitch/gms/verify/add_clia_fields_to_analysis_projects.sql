-- Verify add_clia_fields_to_analysis_projects

BEGIN;

  SELECT is_cle,run_as FROM config.analysis_project WHERE FALSE;

ROLLBACK;
