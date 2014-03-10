-- Revert config_analysis_project

BEGIN;

DROP TABLE IF EXISTS config.analysis_project;

COMMIT;
