-- Revert config_analysis_project_model_bridge

BEGIN;

DROP TABLE IF EXISTS config.analysis_project_model_bridge;

COMMIT;
