-- Verify config_analysis_menu_item_description

BEGIN;

ALTER TABLE config.analysis_menu_item DROP COLUMN description;

ROLLBACK;
