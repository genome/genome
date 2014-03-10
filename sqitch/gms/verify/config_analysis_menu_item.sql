-- Verify config_analysis_menu_item

BEGIN;

SELECT id, created_at, updated_at, name, configuration_set_id
FROM config.analysis_menu_item
WHERE FALSE;

ROLLBACK;
