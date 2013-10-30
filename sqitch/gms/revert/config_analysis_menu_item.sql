-- Revert config_analysis_menu_item

BEGIN;

DROP TABLE IF EXISTS config.analysis_menu_item;

COMMIT;
