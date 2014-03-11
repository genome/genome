-- Verify config_analysis_project_model_bridge

BEGIN;

SELECT id, created_by, created_at, updated_at, analysis_project_id, model_id
FROM config.analysis_project_model_bridge
WHERE FALSE;

ROLLBACK;
