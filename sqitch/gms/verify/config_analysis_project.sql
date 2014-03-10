-- Verify config_analysis_project

BEGIN;

SELECT id, configuration_set_id, created_by, analysis_menu_item_id,
        created_at, updated_at, name, status, model_group_id
FROM config.analysis_project
WHERE FALSE;

ROLLBACK;
