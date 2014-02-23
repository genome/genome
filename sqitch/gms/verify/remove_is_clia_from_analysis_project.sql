-- Verify remove_is_clia_from_analysis_project

BEGIN;

ALTER TABLE config.analysis_project ADD COLUMN is_clia BOOLEAN;

ROLLBACK;
