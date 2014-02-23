-- Verify rename_clia_to_cle_for_analysis_projects

BEGIN;

SELECT is_cle from config.analysis_project WHERE FALSE;

ROLLBACK;
