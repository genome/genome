-- Revert config_analysis_menu_item_index

BEGIN;

DROP INDEX config.analysis_menu_config_set_index;

COMMIT;
